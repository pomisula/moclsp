#include "instance_loader.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using SectionLines = std::unordered_map<std::string, std::vector<std::string>>;

std::vector<double> parse_numbers(std::string line) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line);
    std::vector<double> values;
    double v;
    while (ss >> v) {
        values.push_back(v);
    }
    return values;
}

SectionLines read_sections(const std::filesystem::path& file_path) {
    std::ifstream in(file_path);
    if (!in) {
        throw std::runtime_error("Cannot open instance file: " + file_path.string());
    }

    SectionLines sections;
    std::string current;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!line.empty() && line[0] == '#') {
            continue;
        }
        if (line.front() == '[' && line.back() == ']') {
            current = line.substr(1, line.size() - 2);
            continue;
        }
        if (!current.empty()) {
            sections[current].push_back(line);
        }
    }
    return sections;
}

}  // namespace

Model load_instance_from_flat(const std::filesystem::path& file_path) {
    auto sections = read_sections(file_path);
    Model inst;

    auto index_it = sections.find("INDEX");
    if (index_it == sections.end() || index_it->second.empty()) {
        throw std::runtime_error("INDEX section missing in " + file_path.string());
    }
    auto idx_vals = parse_numbers(index_it->second.front());
    if (idx_vals.size() != 3) {
        throw std::runtime_error("INDEX must have 3 numbers (J, T, M)");
    }
    inst.J = static_cast<int>(idx_vals[0]);
    inst.T = static_cast<int>(idx_vals[1]);
    inst.M = static_cast<int>(idx_vals[2]);

    auto grab_vector = [&](const std::string& key) -> std::vector<double> {
        auto it = sections.find(key);
        if (it == sections.end() || it->second.empty()) {
            return {};
        }
        return parse_numbers(it->second.front());
    };

    inst.holding_cost = grab_vector("LAGKOST");
    inst.initial_inventory = grab_vector("L0");

    if (sections.count("KAPAZ")) {
        for (const auto& ln : sections["KAPAZ"]) {
            inst.capacity.push_back(parse_numbers(ln));
        }
    }

    if (sections.count("P-BEDARF")) {
        for (const auto& ln : sections["P-BEDARF"]) {
            inst.demand.push_back(parse_numbers(ln));
        }
    }

    if (sections.count("RUESTK")) {
        for (const auto& ln : sections["RUESTK"]) {
            auto vals = parse_numbers(ln);
            if (!vals.empty()) {
                inst.setup_cost.push_back(vals.front());
            }
        }
    }

    inst.prod_coeff.assign(inst.M, std::vector<double>(inst.J, 0.0));
    if (sections.count("PRODKOEF")) {
        for (const auto& ln : sections["PRODKOEF"]) {
            auto vals = parse_numbers(ln);
            if (vals.size() != 3) {
                throw std::runtime_error("PRODKOEF lines must have 3 numbers");
            }
            int m = static_cast<int>(vals[0]) - 1;
            int j = static_cast<int>(vals[1]) - 1;
            double a = vals[2];
            if (m < 0 || m >= inst.M || j < 0 || j >= inst.J) {
                throw std::runtime_error("PRODKOEF index out of bounds");
            }
            inst.prod_coeff[m][j] = a;
        }
    }

    inst.setup_time.assign(inst.M, std::vector<double>(inst.J, 0.0));
    if (sections.count("RUESTZ")) {
        for (const auto& ln : sections["RUESTZ"]) {
            auto vals = parse_numbers(ln);
            if (vals.size() != 3) {
                throw std::runtime_error("RUESTZ lines must have 3 numbers");
            }
            int m = static_cast<int>(vals[0]) - 1;
            int j = static_cast<int>(vals[1]) - 1;
            double st = vals[2];
            if (m < 0 || m >= inst.M || j < 0 || j >= inst.J) {
                throw std::runtime_error("RUESTZ index out of bounds");
            }
            inst.setup_time[m][j] = st;
        }
    }

    if (sections.count("UEBER-KS")) {
        for (const auto& ln : sections["UEBER-KS"]) {
            auto vals = parse_numbers(ln);
            if (!vals.empty()) {
                inst.overtime_cost.push_back(vals.front());
            }
        }
    }

    auto expect = [&](int size, int expected, const std::string& label) {
        if (expected != 0 && size != expected) {
            throw std::runtime_error(label + " size mismatch");
        }
    };

    expect(static_cast<int>(inst.holding_cost.size()), inst.J, "LAGKOST");
    expect(static_cast<int>(inst.initial_inventory.size()), inst.J, "L0");
    expect(static_cast<int>(inst.demand.size()), inst.J, "P-BEDARF rows");
    expect(static_cast<int>(inst.capacity.size()), inst.M, "KAPAZ rows");

    for (int j = 0; j < inst.J; ++j) {
        expect(static_cast<int>(inst.demand[j].size()), inst.T, "P-BEDARF cols");
    }
    for (int m = 0; m < inst.M; ++m) {
        expect(static_cast<int>(inst.capacity[m].size()), inst.T, "KAPAZ cols");
    }
    if (!inst.overtime_cost.empty()) {
        expect(static_cast<int>(inst.overtime_cost.size()), inst.M, "UEBER-KS");
    }

    return inst;
}
