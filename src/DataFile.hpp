#include <fstream>
#include <string>

#ifndef DATA_FILE_HPP
#define DATA_FILE_HPP

class DataFile {
public:
  static DataFile *DATA;

  std::ofstream file;

  bool write_file;
  
  DataFile();
  void init();
  void close();

  private : 
  bool open;
};

#endif // DATA_FILE_HPP
