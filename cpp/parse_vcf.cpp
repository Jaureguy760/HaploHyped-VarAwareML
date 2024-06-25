#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "vcfpp.h"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <string>
#include <optional>
#include <memory>
#include <execution>

namespace py = pybind11;
using namespace vcfpp;
using namespace std;

/**
 * Class to handle VCF file operations.
 */
class VCFLoader {
public:
    /**
     * Load VCF data with phased information for a specific sample and chromosome.
     *
     * @param in_vcf Path to the VCF file.
     * @param sample Sample identifier.
     * @param chrom Chromosome identifier.
     * @return A vector of tuples containing SNP data.
     */
    vector<tuple<string, uint32_t, uint32_t, string, string, int8_t, int8_t>> load_vcf(
        const string& in_vcf, const string& sample, const string& chrom) {

        vector<tuple<string, uint32_t, uint32_t, string, string, int8_t, int8_t>> snp_list;
        snp_list.reserve(100000);  // Pre-allocate more memory to avoid frequent reallocations

        try {
            BcfReader vcf(in_vcf, chrom, sample);
            BcfRecord var(vcf.header);
            vector<int> gt(2);  // Assuming diploid, allocate space for 2 genotypes

            while (vcf.getNextVariant(var)) {
                if (!var.isSNP()) {
                    continue;
                }
                var.getGenotypes(gt);
                assert(var.ploidy() == 2);

                string ref = var.REF();
                string alt = var.ALT();

                int8_t phase1 = static_cast<int8_t>(gt[0]);
                int8_t phase2 = static_cast<int8_t>(gt[1]);
                snp_list.emplace_back(
                    var.CHROM(),
                    static_cast<uint32_t>(var.Start()),
                    static_cast<uint32_t>(var.End()),
                    move(ref),
                    move(alt),
                    phase1,
                    phase2
                );
            }
        } catch (const exception& e) {
            cerr << "Error parsing VCF file: " << e.what() << endl;
            throw runtime_error("Error parsing VCF file: " + string(e.what()));
        }

        // Reduced logging: Output only once for total SNPs loaded
        cout << "Loaded " << snp_list.size() << " SNPs for sample " << sample << " and chromosome " << chrom << endl;
        return snp_list;
    }

    /**
     * Load VCF data without sample information for a specific chromosome.
     *
     * @param in_vcf Path to the VCF file.
     * @param chrom Chromosome identifier.
     * @return A vector of tuples containing SNP data.
     */
    vector<tuple<string, uint32_t, uint32_t, string, string>> load_vcf_without_sample(
        const string& in_vcf, const string& chrom) {

        vector<tuple<string, uint32_t, uint32_t, string, string>> snp_list;
        snp_list.reserve(100000);  // Pre-allocate more memory to avoid frequent reallocations

        try {
            BcfReader vcf(in_vcf, chrom);
            BcfRecord var(vcf.header);

            while (vcf.getNextVariant(var)) {
                if (!var.isSNP()) {
                    continue;
                }

                string ref = var.REF();
                string alt = var.ALT();
                snp_list.emplace_back(
                    var.CHROM(),
                    static_cast<uint32_t>(var.Start()),
                    static_cast<uint32_t>(var.End()),
                    move(ref),
                    move(alt)
                );
            }
        } catch (const exception& e) {
            cerr << "Error parsing VCF file: " << e.what() << endl;
            throw runtime_error("Error parsing VCF file: " + string(e.what()));
        }

        // Reduced logging: Output only once for total SNPs loaded
        cout << "Loaded " << snp_list.size() << " SNPs for chromosome " << chrom << endl;
        return snp_list;
    }
};

PYBIND11_MODULE(parse_vcf, m) {
    m.doc() = "Module for parsing VCF files using VCFLoader class";
    py::class_<VCFLoader>(m, "VCFLoader")
        .def(py::init<>())
        .def("load_vcf", &VCFLoader::load_vcf, "Load VCF data with phased information and sample",
             py::arg("in_vcf"), py::arg("sample"), py::arg("chrom") = "")
        .def("load_vcf_without_sample", &VCFLoader::load_vcf_without_sample, "Load VCF data without sample information",
             py::arg("in_vcf"), py::arg("chrom") = "");
}
