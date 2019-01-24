/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "main.hpp"

using namespace std;

int main(int argc, char **argv) {
  if (argc == 2) {
    std::string infilename = argv[1];
    datahandeler(infilename, "out.root");
    // datahandeler2(infilename);
  } else if (argc == 3) {
    std::string infilename = argv[1];
    std::string outfilename = argv[2];
    // std::cerr << RED << "Running SinglePi: \n";
    // SinglePi(infilename, outfilename);
    datahandeler(infilename, outfilename);
    // datahandeler2(infilename);
    std::cerr << RESET << std::endl;
  } else {
    std::cerr << RED << "Error: \n";
    std::cerr << BOLDRED << "\tNeed input file and output file\n";
    std::cerr << RESET << "Usage:\n\t";
    std::cerr << BOLDWHITE << argv[0] << " infile.root outfile.root\n\n";
    std::cerr << RESET << std::endl;
  }

  return 0;
}
