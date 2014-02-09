#include <iostream>
#include "gtest/gtest.h"
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Read.h"

using namespace std;
vector<string> lines;
class VPTest : public ::testing::Test {
protected:
    VPTest() {};
    virtual ~VPTest() {};
    virtual void SetUp() {};
    virtual void TearDown() {};
};

/*
TEST_F(VPTest, BasicParser) {
    Read r1 = Read(lines[0]);
    cout << r1 << endl;
    Read r2 = Read(lines[1]);
    cout << r2 << endl;
    Read r3 = Read(lines[2]);
    cout << r3 << endl;
    Read r4 = Read(lines[3]);
    cout << r4 << endl;
    Read r5 = Read(lines[4]);
    cout << r5 << endl;
}
*/

void should_include_variant(const vector<VariantFromSAM> &variants, VariantFromSAM v) {
    ASSERT_TRUE(std::find(variants.begin(), variants.end(), v) != variants.end())
    << "the variant was not found: " << endl
    << v << endl;
}

TEST_F(VPTest, BASIC) {
    Read r = Read(lines[0]);
    cout << r << endl;
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "A", "T", 1824, 1825));
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "C", "A", 1846, 1847));
}


TEST_F(VPTest, DEL) {
    Read r = Read(lines[1]);
    cout << r << endl;
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::DEL, "GAC", "-", 1826, 1829));
}

TEST_F(VPTest, INS) {
    Read r = Read(lines[2]);
    cout << r << endl;
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::INS, "-", "GTC", 1328, 1328));
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "T", "G", 1338, 1339));
}

TEST_F(VPTest, SOFTCLIP) {
    Read r = Read(lines[3]);
    cout << r << endl;
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::SOFTCLIP, "GCAACAATGG", "TTTTTTTTTT", 1368, 1378));
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "C", "A", 1419, 1420));
}

TEST_F(VPTest, INS2) {
    Read r = Read(lines[4]);
    cout << r << endl;
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::INS, "-", "GTC", 1328, 1328));
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "C", "A", 1304, 1305));
}

TEST_F(VPTest, DEL2) {
    Read r = Read(lines[5]);
    cout << r << endl;
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::DEL, "GAC", "-", 1826, 1829));
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "T", "G", 1780, 1781));
    should_include_variant(r.variants, VariantFromSAM(VariantFromSAM::MISMATCH, "C", "A", 1846, 1847));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        cerr << "please specify the path of data folder" << endl;
        exit(1);
    }
    string test_sam_path = string(argv[1]);
    ifstream varfile((test_sam_path).c_str());
    string line;
    while(getline(varfile, line)) {
        lines.push_back(line);
    }
    return RUN_ALL_TESTS();
}

