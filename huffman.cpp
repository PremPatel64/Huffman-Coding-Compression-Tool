// huffman.cpp
// Build: g++ -std=c++17 -O2 -o huffman huffman.cpp
// Usage:
//   Compress:   ./huffman -c input.txt compressed.huf
//   Decompress: ./huffman -d compressed.huf output.txt

#include <bits/stdc++.h>
using namespace std;

class HuffmanNode {
public:
    uint8_t symbol; // valid only for leaves
    uint64_t freq;
    shared_ptr<HuffmanNode> left;
    shared_ptr<HuffmanNode> right;

    // leaf constructor
    HuffmanNode(uint8_t s, uint64_t f) : symbol(s), freq(f), left(nullptr), right(nullptr) {}
    // internal node constructor
    HuffmanNode(shared_ptr<HuffmanNode> l, shared_ptr<HuffmanNode> r)
        : symbol(0), freq(l->freq + r->freq), left(l), right(r) {}

    bool is_leaf() const { return !left && !right; }
};

struct NodeCompare {
    bool operator()(const shared_ptr<HuffmanNode>& a, const shared_ptr<HuffmanNode>& b) const {
        // priority_queue default is max-heap; we need min-heap based on freq
        return a->freq > b->freq;
    }
};

class HuffmanTree {
public:
    shared_ptr<HuffmanNode> root;

    HuffmanTree() : root(nullptr) {}

    // Build Huffman tree from frequency table (256 entries)
    void build_from_frequencies(const array<uint64_t, 256>& freq) {
        priority_queue<shared_ptr<HuffmanNode>, vector<shared_ptr<HuffmanNode>>, NodeCompare> pq;
        for (size_t i = 0; i < freq.size(); ++i) {
            if (freq[i] > 0) {
                pq.push(make_shared<HuffmanNode>(static_cast<uint8_t>(i), freq[i]));
            }
        }

        if (pq.empty()) {
            root = nullptr;
            return;
        }

        // If only one unique symbol, duplicate it to create a two-leaf tree
        if (pq.size() == 1) {
            auto only = pq.top(); pq.pop();
            auto dup = make_shared<HuffmanNode>(only->symbol, only->freq);
            root = make_shared<HuffmanNode>(only, dup);
            return;
        }

        while (pq.size() > 1) {
            auto l = pq.top(); pq.pop();
            auto r = pq.top(); pq.pop();
            pq.push(make_shared<HuffmanNode>(l, r));
        }
        root = pq.top();
    }

    // Generate binary codes ("0"/"1" strings) for each byte value present
    void build_codes(unordered_map<uint8_t, string>& out_codes) const {
        out_codes.clear();
        if (!root) return;
        string cur;
        traverse(root, cur, out_codes);
    }

private:
    void traverse(shared_ptr<HuffmanNode> node, string& cur, unordered_map<uint8_t, string>& out_codes) const {
        if (!node) return;
        if (node->is_leaf()) {
            out_codes[node->symbol] = cur.empty() ? "0" : cur; // single-symbol case -> "0"
            return;
        }
        cur.push_back('0');
        traverse(node->left, cur, out_codes);
        cur.pop_back();

        cur.push_back('1');
        traverse(node->right, cur, out_codes);
        cur.pop_back();
    }
};

class HuffmanCodec {
public:
    // Compress file: returns pair{original_size, compressed_file_size}, and sets elapsed_seconds
    static pair<uint64_t, uint64_t> compress_file(const string& in_path, const string& out_path, double& elapsed_seconds) {
        auto t0 = chrono::high_resolution_clock::now();

        // Read entire input file into memory
        vector<uint8_t> input_bytes = read_all_bytes(in_path);
        uint64_t original_size = input_bytes.size();

        // Frequency table
        array<uint64_t, 256> freq{};
        freq.fill(0);
        for (uint8_t b : input_bytes) ++freq[b];

        // Build Huffman tree and codes
        HuffmanTree tree;
        tree.build_from_frequencies(freq);
        unordered_map<uint8_t, string> codes;
        tree.build_codes(codes);

        // Encode bitstream into bytes (MSB-first)
        vector<uint8_t> compressed_bytes;
        uint8_t cur_byte = 0;
        int bit_count = 0; // how many bits filled in cur_byte (0..7)
        for (uint8_t b : input_bytes) {
            const string& code = codes[b];
            for (char ch : code) {
                cur_byte <<= 1;
                if (ch == '1') cur_byte |= 1;
                ++bit_count;
                if (bit_count == 8) {
                    compressed_bytes.push_back(cur_byte);
                    cur_byte = 0;
                    bit_count = 0;
                }
            }
        }
        uint8_t padding = 0;
        if (bit_count != 0) {
            // shift remaining bits to MSB side and pad LSB with zeros
            cur_byte <<= (8 - bit_count);
            compressed_bytes.push_back(cur_byte);
            padding = static_cast<uint8_t>(8 - bit_count);
        }

        // Write header + payload
        ofstream out(out_path, ios::binary);
        if (!out) throw runtime_error("Failed to open output file: " + out_path);

        // Header: magic (4 bytes) + original_size + unique_count + (symbol,freq)* + padding
        const char magic[4] = { 'H', 'U', 'F', '1' };
        out.write(magic, 4);
        write_uint64(out, original_size);

        uint16_t unique_count = 0;
        for (auto f : freq) if (f > 0) ++unique_count;
        write_uint16(out, unique_count);

        for (size_t i = 0; i < freq.size(); ++i) {
            if (freq[i] > 0) {
                uint8_t sym = static_cast<uint8_t>(i);
                out.write(reinterpret_cast<const char*>(&sym), sizeof(sym));
                write_uint64(out, freq[i]);
            }
        }

        out.write(reinterpret_cast<const char*>(&padding), sizeof(padding));
        if (!compressed_bytes.empty()) {
            out.write(reinterpret_cast<const char*>(compressed_bytes.data()), compressed_bytes.size());
        }
        out.close();

        auto t1 = chrono::high_resolution_clock::now();
        elapsed_seconds = chrono::duration<double>(t1 - t0).count();

        uint64_t compressed_size = file_size(out_path);
        return { original_size, compressed_size };
    }

    // Decompress file; returns original_size (as read from header) and sets elapsed_seconds
    static uint64_t decompress_file(const string& in_path, const string& out_path, double& elapsed_seconds) {
        auto t0 = chrono::high_resolution_clock::now();

        ifstream in(in_path, ios::binary);
        if (!in) throw runtime_error("Failed to open input file: " + in_path);

        // Read and validate magic
        char magic[4];
        in.read(magic, 4);
        if (in.gcount() != 4 || magic[0] != 'H' || magic[1] != 'U' || magic[2] != 'F' || magic[3] != '1') {
            throw runtime_error("Invalid or unsupported compressed file format (bad magic).");
        }

        uint64_t original_size = read_uint64(in);
        uint16_t unique_count = read_uint16(in);

        array<uint64_t, 256> freq{};
        freq.fill(0);
        for (uint16_t i = 0; i < unique_count; ++i) {
            uint8_t sym = 0;
            in.read(reinterpret_cast<char*>(&sym), sizeof(sym));
            uint64_t f = read_uint64(in);
            freq[sym] = f;
        }

        uint8_t padding = 0;
        in.read(reinterpret_cast<char*>(&padding), sizeof(padding));

        // Read remaining bytes into compressed_bytes
        vector<uint8_t> compressed_bytes;
        {
            // Move to current pos and read rest
            istreambuf_iterator<char> it(in), end;
            for (; it != end; ++it) compressed_bytes.push_back(static_cast<uint8_t>(*it));
        }
        in.close();

        // Rebuild tree
        HuffmanTree tree;
        tree.build_from_frequencies(freq);
        if (!tree.root) {
            // empty original file
            ofstream out(out_path, ios::binary);
            out.close();
            elapsed_seconds = chrono::duration<double>(chrono::high_resolution_clock::now() - t0).count();
            return 0;
        }

        // Decode bits by traversing tree
        vector<uint8_t> output_bytes;
        output_bytes.reserve(static_cast<size_t>(original_size));

        size_t total_bits = compressed_bytes.size() * 8;
        if (!compressed_bytes.empty()) total_bits -= padding;

        size_t bit_index = 0;
        shared_ptr<HuffmanNode> node = tree.root;

        for (size_t byte_idx = 0; byte_idx < compressed_bytes.size() && output_bytes.size() < original_size; ++byte_idx) {
            uint8_t cur = compressed_bytes[byte_idx];
            for (int bitpos = 7; bitpos >= 0; --bitpos) {
                if (bit_index >= total_bits) break;
                bool bit = ((cur >> bitpos) & 1) != 0;
                node = bit ? node->right : node->left;
                if (!node) throw runtime_error("Corrupt compressed data (unexpected null node).");
                if (node->is_leaf()) {
                    output_bytes.push_back(node->symbol);
                    node = tree.root;
                    if (output_bytes.size() >= original_size) break;
                }
                ++bit_index;
            }
        }

        // Write decompressed bytes to output file
        ofstream out(out_path, ios::binary);
        if (!out) throw runtime_error("Failed to open output file for writing: " + out_path);
        if (!output_bytes.empty()) out.write(reinterpret_cast<const char*>(output_bytes.data()), output_bytes.size());
        out.close();

        auto t1 = chrono::high_resolution_clock::now();
        elapsed_seconds = chrono::duration<double>(t1 - t0).count();
        return original_size;
    }

private:
    // ----------------- Helper I/O and utils -----------------
    static vector<uint8_t> read_all_bytes(const string& path) {
        ifstream in(path, ios::binary | ios::ate);
        if (!in) throw runtime_error("Failed to open input file: " + path);
        streamsize size = in.tellg();
        if (size < 0) size = 0;
        in.seekg(0, ios::beg);
        vector<uint8_t> buf(static_cast<size_t>(size));
        if (size > 0) in.read(reinterpret_cast<char*>(buf.data()), size);
        return buf;
    }

    static uint64_t file_size(const string& path) {
        ifstream in(path, ios::binary | ios::ate);
        if (!in) return 0;
        return static_cast<uint64_t>(in.tellg());
    }

    static void write_uint64(ofstream& out, uint64_t v) {
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
    static uint64_t read_uint64(ifstream& in) {
        uint64_t v = 0;
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        return v;
    }
    static void write_uint16(ofstream& out, uint16_t v) {
        out.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
    static uint16_t read_uint16(ifstream& in) {
        uint16_t v = 0;
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        return v;
    }
};

static void print_usage(const char* prog) {
    cerr << "Usage:\n"
         << "  Compress:   " << prog << " -c input.txt output.huf\n"
         << "  Decompress: " << prog << " -d input.huf output.txt\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        if (argc != 4) {
            print_usage(argv[0]);
            return 1;
        }

        string mode = argv[1];
        string in_path = argv[2];
        string out_path = argv[3];

        if (mode == "-c") {
            double elapsed = 0.0;
            auto sizes = HuffmanCodec::compress_file(in_path, out_path, elapsed);
            uint64_t original = sizes.first;
            uint64_t compressed = sizes.second;

            cout.setf(ios::fixed);
            cout << setprecision(6);
            cout << "Compression finished.\n";
            cout << " Original size:   " << original << " bytes\n";
            cout << " Compressed size: " << compressed << " bytes\n";
            if (original == 0) {
                cout << " Compression ratio: N/A (empty input)\n";
            } else {
                double ratio = static_cast<double>(compressed) / static_cast<double>(original);
                double savings = (1.0 - ratio) * 100.0;
                cout << " Compression ratio: " << ratio << " (compressed / original)\n";
                cout << " Space savings:      " << savings << " %\n";
            }
            cout << " Time elapsed: " << elapsed << " seconds\n";
        } else if (mode == "-d") {
            double elapsed = 0.0;
            uint64_t out_size = HuffmanCodec::decompress_file(in_path, out_path, elapsed);
            cout << "Decompression finished.\n";
            cout << " Output size: " << out_size << " bytes\n";
            cout << " Time elapsed: " << elapsed << " seconds\n";
        } else {
            print_usage(argv[0]);
            return 1;
        }
    } catch (const exception& ex) {
        cerr << "ERROR: " << ex.what() << "\n";
        return 2;
    }

    return 0;
}
