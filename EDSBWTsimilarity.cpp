//
//  EDSBWTsimilarity
//
//  Programma unico per calcolare la similarita simmetrica tra due EDS:
//  1. costruisce o riusa l indice EDS-BWT della prima EDS;
//  2. costruisce o riusa l indice EDS-BWT della seconda EDS;
//  3. calcola A(P1,P2), cercando P1 nell indice di P2;
//  4. calcola A(P2,P1), cercando P2 nell indice di P1;
//  5. restituisce la media aritmetica.
//

#include <cerrno>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include "EDSBWTsearch.hpp"
#include "Parameters.h"

using namespace std;

static bool fileExists(const string& path)
{
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

static bool fileNewerThan(const string& file, const string& reference)
{
    struct stat fileInfo;
    struct stat referenceInfo;

    if (stat(file.c_str(), &fileInfo) != 0 || stat(reference.c_str(), &referenceInfo) != 0) {
        return false;
    }

    return fileInfo.st_mtime >= referenceInfo.st_mtime;
}

static bool ensureDirectory(const string& path)
{
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        return S_ISDIR(info.st_mode);
    }

    if (mkdir(path.c_str(), 0775) != 0) {
        cerr << "Could not create directory " << path << ": " << strerror(errno) << endl;
        return false;
    }

    return true;
}

static string basenameWithoutEdsExtension(const string& path)
{
    size_t slash = path.find_last_of("/");
    string name = slash == string::npos ? path : path.substr(slash + 1);

    const string ext = ".eds";
    if (name.size() >= ext.size() && name.substr(name.size() - ext.size()) == ext) {
        return name.substr(0, name.size() - ext.size());
    }

    return name;
}

static string indexBaseForEDS(const string& edsPath)
{
    return string("edsbwt_form/") + basenameWithoutEdsExtension(edsPath);
}

static int runCommand(const vector<string>& args)
{
    if (args.empty()) {
        cerr << "Internal error: empty command" << endl;
        return 1;
    }

    vector<char*> argv;
    for (size_t i = 0; i < args.size(); i++) {
        argv.push_back(const_cast<char*>(args[i].c_str()));
    }
    argv.push_back(NULL);

    pid_t pid = fork();
    if (pid < 0) {
        cerr << "fork failed while running " << args[0] << ": " << strerror(errno) << endl;
        return 1;
    }

    if (pid == 0) {
        execvp(argv[0], argv.data());
        cerr << "exec failed for " << args[0] << ": " << strerror(errno) << endl;
        _exit(127);
    }

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        cerr << "waitpid failed for " << args[0] << ": " << strerror(errno) << endl;
        return 1;
    }

    if (WIFEXITED(status)) {
        return WEXITSTATUS(status);
    }

    if (WIFSIGNALED(status)) {
        cerr << args[0] << " terminated by signal " << WTERMSIG(status) << endl;
        return 128 + WTERMSIG(status);
    }

    return 1;
}

static void removeIfExists(const string& path)
{
    if (fileExists(path) && remove(path.c_str()) != 0) {
        cerr << "Warning: could not remove " << path << ": " << strerror(errno) << endl;
    }
}

static void cleanupIndexFiles(const string& outputBase)
{
    removeIfExists(outputBase + ".fasta");
    removeIfExists(outputBase + ".len");
    removeIfExists(outputBase + ".info");
    removeIfExists(outputBase + ".empty.info");
    removeIfExists(outputBase + ".bwt");
    removeIfExists(outputBase + ".4.da");
    removeIfExists(outputBase + ".ebwt");
    removeIfExists(outputBase + ".lcp");
    removeIfExists(outputBase + ".da");
    removeIfExists(outputBase + ".posSA");
    removeIfExists(outputBase + ".SAP");
    removeIfExists(outputBase + ".bitvector");
    removeIfExists(outputBase + "_info.aux");
    removeIfExists(outputBase + "_alpha.txt");
    removeIfExists(outputBase + "_tableOcc.txt");

    for (int i = 0; i < 256; i++) {
        string suffix = to_string(i) + ".aux";
        removeIfExists(outputBase + "_bwt_" + suffix);
        removeIfExists(outputBase + "_bv_" + suffix);
    }
}

static bool indexIsReady(const string& edsPath, const string& outputBase)
{
    return fileExists(outputBase + ".ebwt") &&
           fileExists(outputBase + "_info.aux") &&
           fileNewerThan(outputBase + ".ebwt", edsPath) &&
           fileNewerThan(outputBase + "_info.aux", edsPath);
}

static int buildEDSBWTIndex(const string& edsPath, const string& outputBase)
{
    cout << "\n=== EDS-BWT index ===" << endl;
    cout << "EDS input: " << edsPath << endl;
    cout << "Index base: " << outputBase << endl;

    if (!fileExists(edsPath)) {
        cerr << "Input EDS file does not exist: " << edsPath << endl;
        return 1;
    }

    if (!ensureDirectory("edsbwt_form")) {
        return 1;
    }

    if (indexIsReady(edsPath, outputBase)) {
        cout << "Index already available, reusing it." << endl;
        return 0;
    }

    cout << "Index missing or older than EDS, rebuilding it." << endl;
    cleanupIndexFiles(outputBase);
    int status = runCommand(vector<string>{"./eds_to_fasta", edsPath, outputBase});
    if (!fileExists(outputBase + ".fasta")) {
        cerr << "eds_to_fasta failed for " << edsPath << " with exit code " << status << endl;
        cleanupIndexFiles(outputBase);
        return status == 0 ? 1 : status;
    }

    status = runCommand(vector<string>{"gsufsort/gsufsort", outputBase + ".fasta",
                                       "--da", "--bwt", "--output", outputBase});
    if (!fileExists(outputBase + ".bwt") || !fileExists(outputBase + ".4.da")) {
        cerr << "gsufsort failed for " << outputBase << " with exit code " << status << endl;
        cleanupIndexFiles(outputBase);
        return status == 0 ? 1 : status;
    }

    removeIfExists(outputBase + ".fasta");
    removeIfExists(outputBase + ".len");
    removeIfExists(outputBase + ".info");

    status = runCommand(vector<string>{"./da_to_everything", outputBase});
    removeIfExists(outputBase + ".empty.info");
    removeIfExists(outputBase + ".bwt");
    removeIfExists(outputBase + ".4.da");

    if (!fileExists(outputBase + "_info.aux")) {
        cerr << "da_to_everything failed for " << outputBase << " with exit code " << status << endl;
        cleanupIndexFiles(outputBase);
        return status == 0 ? 1 : status;
    }

    return 0;
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " P1.eds P2.eds" << endl;
        cerr << "where P1.eds and P2.eds are two elastic-degenerate strings" << endl;
        return 1;
    }

    string p1EDS = argv[1];
    string p2EDS = argv[2];
    string p1Index = indexBaseForEDS(p1EDS);
    string p2Index = indexBaseForEDS(p2EDS);

    int status = buildEDSBWTIndex(p1EDS, p1Index);
    if (status != 0) {
        return status;
    }

    status = buildEDSBWTIndex(p2EDS, p2Index);
    if (status != 0) {
        return status;
    }

    cout << "\n=== Computing directional similarities ===" << endl;

    float a12 = 0.0f;
    float a21 = 0.0f;

    {
        EDSBWT indexP1(p1Index, MODE, 1);
        EDSBWT indexP2(p2Index, MODE, 1);

        cout << "\n--- A(P1,P2): searching P1 in index(P2) ---" << endl;
        a12 = indexP2.computeSimilarityFromEDS(p1EDS);

        cout << "\n--- A(P2,P1): searching P2 in index(P1) ---" << endl;
        a21 = indexP1.computeSimilarityFromEDS(p2EDS);
    }

    float finalSimilarity = (a12 + a21) / 2.0f;

    cout << "\n=== Final similarity ===" << endl;
    cout << "A(P1,P2): " << a12 << endl;
    cout << "A(P2,P1): " << a21 << endl;
    cout << "Similarity: " << finalSimilarity << endl;

    cout << "\nIndexes kept in edsbwt_form/." << endl;
    return 0;
}
