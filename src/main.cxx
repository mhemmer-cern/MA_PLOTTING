#include "../include/main.hpp"

int main(int argc, char const *argv[]) {

  if(argc == 10){
    plotting(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], atoi(argv[9]));
  }
  else{
    std::cout << std::string(80, '_') << std::endl;
    std::cout << "|Number of Arguments given: " << argc << std::endl;
    std::cout << "| WRONG NUMBER OF ARGUMENTS! EXIT CODE 3 " << std::endl;
    std::cout << std::string(80, '_') << std::endl;
    exit(3);
  }

  return 0;
}
