/*
 * ConfigFile.h
 *   Thanks: this file was Obtained from http://www.adp-gmbh.ch/cpp/
 *   It is used to read a property file !!
 */

#ifndef CONFIGFILE_H_
#define CONFIGFILE_H_

#include <string>
#include <map>

#include "Chameleon.h"

class ConfigFile {
  std::map<std::string,Chameleon> content_;

public:
  ConfigFile(std::string const& configFile);

  Chameleon const& Value(std::string const& section, std::string const& entry) const;

  Chameleon const& Value(std::string const& section, std::string const& entry, double value);
  Chameleon const& Value(std::string const& section, std::string const& entry, std::string const& value);
};

#endif /* CONFIGFILE_H_ */
