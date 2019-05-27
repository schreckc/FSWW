#include "DataFile.hpp"

DataFile *DataFile::DATA(new DataFile());

DataFile::DataFile() {
  open = false;
  write_file = false;
}

void DataFile::init() {
  if (open) {
    close();
  }
  open = true;
  if (write_file) {
    file.open("data/file.dat");
  }
}

void DataFile::close() {
  if (open) {
    open = false;
    if (write_file) {
      file.close();
    }
  }
}
