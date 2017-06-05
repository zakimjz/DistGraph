/*
 * Copyright (C) 2014 Robert Kessl
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "logger.hpp"
#include <time.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <utils.hpp>

using utils::split;

using namespace std;


// static variables
map<string, Logger *> *Logger::loggers = 0;
int Logger::log_message_number = 0;

log_level_t Logger::default_log_level = DEFAULT_LOG_LEVEL;
bool Logger::log_level_read = false;

int Logger::processor_rank_static = 0;
int Logger::processor_number_static = 0;
bool Logger::runs_in_parallel = false;
streambuf *Logger::global_file_buffer = 0;


map<string, log_level_t> *Logger::log_levels = 0;



/** Manipulator that logs header of the log entry
 */
Logger & head(Logger& os)
{
  return os.loghead();
}

/** Manipulator that starts new entry
 */
Logger & ele(Logger & os)
{
  if(os.get_detail_flag() == false) {
    os.put(os.widen('\n')).flush();
  }
  os.too_detailed(false);
  os.writing_le(false);
  return os;
}

void Logger::init()
{
  this->logger_name = "";
  this->loglevel = DEFAULT_LOG_LEVEL;
  this->local_timestamp = 0;
  this->log_head_format = DEFAULT_LOG_FORMAT;
  this->_writing_le = false;
  stream_buffer = 0;
  _too_detailed = false;
  _this_entry_loglevel = LG_INFO;
  have_source_information = false;
  _writing_le = false;
}

string Logger::process_filename(const string &fname)
{
  size_t last_slash = fname.find_last_of('/');
  if(last_slash == string::npos) return fname;
  string result = fname.substr(last_slash + 1);
  return result;
}

Logger::Logger(const Logger &log) : ostream(log.rdbuf())
{
  init();
  logger_name = log.logger_name;
  loglevel = log.loglevel;
  local_timestamp = log.local_timestamp;
  log_head_format = log.log_head_format;
  stream_buffer = log.stream_buffer;
  _too_detailed = log._too_detailed;
  _this_entry_loglevel = log._this_entry_loglevel;

  processor_rank = 0;
  processor_number = 0;
}

Logger::Logger(string logger_name, const char *log_head_format, log_level_t loglevel, int default_log_level, const ostream &output)
  : ostream(output.rdbuf())
{
  init();
  this->logger_name = logger_name;
  this->loglevel = loglevel;
  this->log_head_format = (const char *)log_head_format;
  stream_buffer = output.rdbuf();

  processor_rank = 0;
  processor_number = 0;
} // Logger


Logger::Logger(string logger_name, const char *log_head_format, log_level_t loglevel, int default_log_level, const char* filename)
  : ostream(new filebuf())
{
  init();
  this->logger_name = logger_name;
  this->loglevel = loglevel;
  this->log_head_format = log_head_format;
  ((filebuf*)rdbuf())->open(filename, ios_base::out | ios_base::trunc);

  processor_rank = 0;
  processor_number = 0;
}

Logger::Logger(string logger_name, const char *log_head_format, log_level_t loglevel, int default_log_level, streambuf *obuf)
  : ostream(obuf)
{
  init();
  this->logger_name = logger_name;
  this->loglevel = loglevel;
  this->log_head_format = (const char *)log_head_format;
  stream_buffer = obuf;

  processor_rank = 0;
  processor_number = 0;
}


Logger::~Logger()
{
  if(stream_buffer != global_file_buffer && stream_buffer != cout.rdbuf()) delete stream_buffer;
}

/** supported parameters:
 *    output.format
 *    output.
 *    logger.defaultlevel
 *    component.
 *    logger.name
 *    logger.mark
 */
void Logger::set_parameter(char *pname, char *value)
{
}


/** Parameters in the head log format
 *  %p - processor rank
 *  %P - processor number
 *  %l - log level (text)
 *  %n - logger name
 *  %T - timestamp
 *  %t - number of log message
 *  %a - thread id
 *  %f - filename
 *  %L - line number
 */
Logger & Logger::loghead()
{
  if(_too_detailed == true) return *this;
  string result;
  unsigned int i;
  //  char *begin;
  char pnum[100];

  time_t ctime;
  tm *broken_ctime;

  _writing_le = true;

  for(i = 0; i < log_head_format.size(); i++) {
    if(log_head_format[i] == '\\' && log_head_format[i + 1] == '%') {
      i++;
      result = result + log_head_format[i];
    } else if(log_head_format[i] == '%') {
      i++;
      switch(log_head_format[i]) {
      case 'P':
        sprintf(pnum, "%d", processor_rank);
        result = result + pnum;
        break;
      case 'p':
        sprintf(pnum, "%d", processor_number);
        result = result + pnum;
        break;
      case 't':
        sprintf(pnum, "%d", log_message_number);
        log_message_number++;
        result =  result + pnum;
        break;
      case 'l':
        result = result + this_entry_loglevel();
        break;
      case 'n':
        result = result + logger_name;
        break;
      case 'T':
        time((time_t *)&ctime);
        broken_ctime = localtime((time_t *)&ctime);
        strftime(pnum, 100, "%a %b %d %H:%M:%S", broken_ctime);
        result = result + pnum;
        //strftime(pnum, 100, "%c",
        break;
      /*      case 'a':
         //sprintf(pnum, "%ld", pthread_self());
         result = result + pnum;
         break;
       */
      case 'L':
        if(have_source_information) {
          sprintf(pnum, "%d", line);
          result = result + pnum;
        } else result += "<line>";
        break;
      case 'f':
        if(have_source_information) {
          result = result + process_filename(file);
        } else result = result + "<file>";
        break;
      default:
        break;
      } // switch
    } else {
      result = result + log_head_format[i];
    } //if-else
  } // for i
  (*this) << result.c_str();
  return *this;
} // loghead



Logger *Logger::get_logger(string component, string format)
{
  if(global_file_buffer == 0) {
    return get_logger(component, format, cout.rdbuf());
  } else {
    return get_logger(component, format, global_file_buffer);
  }
}


Logger *Logger::get_logger(string component, string format, string filename)
{

  streambuf *buf = 0;
  if(global_file_buffer == 0) {
    buf = new filebuf();
    ((filebuf*)buf)->open(filename.c_str(), ios_base::out | ios_base::trunc);
  } else {
    buf = global_file_buffer;
  }
  return get_logger(component, format, buf);
}


Logger *Logger::get_logger(string component, string format, streambuf *obuf)
{
  log_level_t ll = get_loglevel_from_env(component);
  if(loggers == 0) loggers = new map<string, Logger *>;

  Logger *logger = 0;

  if(loggers->find(component) == loggers->end()) {
    logger = loggers->insert(pair<string,Logger*>(component, new Logger(component, format.c_str(), ll, 0, obuf))).first->second;
  } else {
    logger = loggers->find(component)->second;
  }

  if(runs_in_parallel) {
    logger->set_processor_rank(processor_rank_static);
    logger->set_processor_number(processor_number_static);
  }

  return logger;
}


Logger *Logger::get_logger(string component)
{
  log_level_t ll = get_loglevel_from_env(component);
  if(loggers == 0) loggers = new map<string, Logger *>;

  if(global_file_buffer != 0) {
    return get_logger(component, DEFAULT_LOG_FORMAT, global_file_buffer);
  } else {
    return get_logger(component, DEFAULT_LOG_FORMAT, cout.rdbuf());
  }
}



void Logger::free_loggers()
{
  for(map<string, Logger *>::iterator itL = loggers->begin(); itL != loggers->end(); itL++) {
    delete itL->second;
  } // for
  delete loggers;
  delete global_file_buffer;
  delete log_levels;
} // free_loggers



log_level_t Logger::get_loglevel_from_env(string logger) {
  if(!log_level_read) {
    char * ll_env = 0;
    ll_env = getenv("LOG_LEVEL");
    if(ll_env == 0) return DEFAULT_LOG_LEVEL;
    std::vector<string> tokens;
    split(ll_env, tokens, ",");

    if(log_levels == 0) log_levels = new map<string, log_level_t>();
    for(std::vector<string>::iterator itV = tokens.begin();
        itV != tokens.end();
        itV++) {
      vector<string> tmp;
      split(*itV, tmp, ":");
      if(tmp.size() == 1) {
        default_log_level = get_log_level(tmp[0].c_str());
        cout << "default log level=" << get_log_level(default_log_level) << endl;
        continue;
      }
      if(tmp.size() != 2) {
        cout << "error while parsing: " << *itV << endl;
        continue;
      }
      cout << "logger: " << tmp[0] << " with level=" << tmp[1] << endl;
      log_level_t ll = get_log_level(tmp[1].c_str());
      cout << "log level: " << ll << endl;
      log_levels->insert(map<string, log_level_t>::value_type(tmp[0], ll));
    }
    log_level_read = true;
  }
  if(log_levels->find(logger) == log_levels->end()) return default_log_level;
  return (*log_levels)[logger];
}


void Logger::disable_loggers_from_env()
{

}



void Logger::set_streambuf(streambuf *buf)
{
  cout << buf << endl;
  ios::rdbuf(buf);
  stream_buffer = buf;
}




void Logger::set_parallel_info(int processor_rank, int processor_number, bool set_to_existing_loggers)
{
  processor_rank_static = processor_rank;
  processor_number_static = processor_number;
  runs_in_parallel = true; //enable_parallel;
  if(set_to_existing_loggers) {
    for(map<string, Logger *>::iterator it = loggers->begin(); it != loggers->end(); it++) {
      it->second->set_processor_rank(processor_rank);
      it->second->set_processor_number(processor_number);
    } // for it
  } // if(set_to_loggers)
}


void Logger::set_parallel_info(int processor_rank, int processor_number, std::string filename_prefix)
{
  set_parallel_info(processor_rank, processor_number);

  std::stringstream ss;
  ss << filename_prefix.c_str() << processor_rank;
  std::filebuf *fbuf = new std::filebuf();
  fbuf->open(ss.str().c_str(), std::ios::out | std::ios::trunc);
  set_global_streambuf(fbuf, true);
}

void Logger::set_global_streambuf(streambuf *buf, bool set_to_loggers)
{
  if(set_to_loggers) {
    for(map<string, Logger *>::iterator it = loggers->begin(); it != loggers->end(); it++) {
      it->second->set_streambuf(buf);
    } // for it
  } // if(set_to_loggers)

  delete global_file_buffer;
  global_file_buffer = buf;
}

