#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <set>
#include <iterator>

// Functions provided by this header are for writing/reading tsv_tables

typedef std::vector<std::vector<std::string>> tsv_table;

inline tsv_table read_tsv(const std::string & fp) {
    tsv_table table;
    std::ifstream file(fp);
    int pos = file.tellg();
    table.reserve(std::count(std::istreambuf_iterator<char>(file),
                             std::istreambuf_iterator<char>(), '\n') + 1);
    file.seekg(pos);
    for (std::string s; std::getline(file, s);) {

        std::istringstream ss(s);
        pos = ss.tellg();
        table.emplace_back(std::vector<std::string>(std::count(std::istreambuf_iterator<char>(ss),
                                                               std::istreambuf_iterator<char>(), '\t') + 1));
        ss.seekg(pos);
        int i = 0;
        for (std::string entry; std::getline(ss, entry, '\t');) {
            table.back()[i] = entry;
            i++;
        }

    }

    uint32_t length = 0;
    for (auto & v : table) if (v.size() > length) length = v.size();
    for (auto& v : table) v.resize(length);
    file.close();
    return table;
}

typedef std::vector<std::string> tsv_table_delineated;



inline tsv_table_delineated read_tsv_delineated(const std::string & fp) {
    tsv_table_delineated table;
    std::ifstream file(fp);
    int pos = file.tellg();
    table.reserve(std::count(std::istreambuf_iterator<char>(file),
                             std::istreambuf_iterator<char>(), '\n') + 1);
    file.seekg(pos);
    for (std::string s; std::getline(file, s);) {
        table.push_back(s);
    }

    file.close();

    return table;
}

inline void write_tsv_table (const tsv_table table, const std::string & fp) {
    std::ofstream file(fp);
    int i = 0;
    for (; i < table.size() - 1; i++) {
        auto it = table[i].begin();
        for (; it != --table[i].end(); ++it) {
            file << *it << '\t';
        }
        file << *it << '\n';
    }
    auto it = table[i].begin();
    for (; it != --table[i].end(); ++it) {
        file << *it << '\t';
    }
    file << *it;
    file.close();
} 

inline std::string join_container(const std::vector<int> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::ostringstream oss;
    std::copy(data.begin(), data.end() - 1, std::ostream_iterator<int>(oss, delimiter));
    std::string result(oss.str());
    result += std::to_string(data.back());
    return result;

}

inline std::string join_container(const std::deque<int> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::ostringstream oss;
    std::copy(data.begin(), data.end() - 1, std::ostream_iterator<int>(oss, delimiter));
    std::string result(oss.str());
    result += std::to_string(data.back());
    return result;

}

inline std::string join_container(const std::vector<uint32_t> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::ostringstream oss;
    std::copy(data.begin(), data.end() - 1, std::ostream_iterator<int>(oss, delimiter));
    std::string result(oss.str());
    result += std::to_string(data.back());
    return result;

}

inline std::string join_container(const std::deque<uint32_t> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::ostringstream oss;
    std::copy(data.begin(), data.end() - 1, std::ostream_iterator<int>(oss, delimiter));
    std::string result(oss.str());
    result += std::to_string(data.back());
    return result;

}

inline std::string join_container(const std::deque<uint8_t> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::ostringstream oss;
    std::copy(data.begin(), data.end() - 1, std::ostream_iterator<int>(oss, delimiter));
    std::string result(oss.str());
    result += std::to_string(data.back());
    return result;

}

inline std::string join_container(const std::vector<std::string> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::string result;

    for (auto it = data.begin(); it != --data.end(); ++it) result += *it + delimiter;

    result += data.back();
    return result;

}

inline std::string join_container(const std::deque<std::string> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::string result;

    if (data.empty()) return {};
    for (auto it = data.begin(); it != --data.end(); ++it) result += *it + delimiter;

    result += data.back();
    return result;

}

inline std::string join_container(const std::set<std::string> & data, const char * delimiter) {

    if (data.empty()) return {};
    std::string result;
    for (auto it = data.begin(); it != --data.end(); ++it) result += *it + delimiter;

    result += *--data.end();
    return result;

}