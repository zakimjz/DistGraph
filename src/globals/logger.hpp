/*
 * Copyright (C) 2014 Robert Kessl
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_
#include <fstream>
#include <map>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdexcept>

#include <iterator>
#include <sstream>
#include <stdlib.h>
#include <string.h>


#define DEFAULT_LOG_FORMAT "%P [%t] [%l] [%T] [%n] [%f@%L] - "
#define DEFAULT_LOG_LEVEL LG_TRACE

typedef enum {LG_OFF, LG_CRITICAL_ERROR, LG_ERROR, LG_WARNING, LG_INFO, LG_DEBUG, LG_TRACE, LG_TRACE1, LG_TRACE2, LG_TRACE3, LG_TRACE4, LG_TRACE5} log_level_t;


using std::string;
using std::streambuf;
using std::filebuf;
using std::map;
using std::ostream;
using std::endl;
using std::cout;
using std::pair;

/* Co by takovy logger mel umet ?
 * vystup by mel automaticky labelovat cislem procesoru, casovym razitkem a cislem threadu
 * mel by podporovat komponenty
 * mel by podporovat streamy i neco jako printf
 * redirect do souboru, zapis do vice souboru najednou (?)
 * nastavovani parametru logu pres funkci set_parameter a info o parametru pres get_parameter
 * log levely
 * thread safe !
 *
 * pozor ja mam:
 * a) loglevel - tedy cislo, ktere mi rika jak detailne mam logovat
 * b) typ logovaci hlasky ! (CRITICAL_ERROR, ERROR, WARNING, INFO, CONFIG, DEBUG1, DEBUG2, ...)
 * CRITICAL_ERROR -
 * ERROR - chyba ze ktere se obcas lze zotavit
 * WARNING - chybne parametry, spatna konfigurace
 * ^--------- nema smysl odlisovat detaily, napr error level 1, error level 2 nebo critical error 1, 2 apod.
 * INFO - informace pro usera, mozna by melo byt vice urovni INFO
 * DEBUG - detailni informace o behu
 * TRACE - vypisuje kazdy krok
 */

/* Jeste poznamky k HashTrie:
 * pokud to udelam jako sablonu kde prvek te usporadane mnoziny,
 * kterou ukladam do HashTrie, bude jako parametr tak asi uz nebudu
 * moci nijak rozume implementovat odlisnost mezi vnitrinmi uzly a
 * listy. Myslim pomoci sablon.  Respektive, bude to hodne slozity.
 *
 */

#define _CRITICAL_ERROR() critical_error(__LINE__, __FILE__)
#define _ERROR() error(__LINE__, __FILE__)
#define _WARNING() warning(__LINE__, __FILE__)
#define _INFO() info(__LINE__, __FILE__)
#define _DEBUG() debug(__LINE__, __FILE__)
#define _TRACE() trace(__LINE__, __FILE__)
#define _TRACE1() trace1(__LINE__, __FILE__)
#define _TRACE2() trace2(__LINE__, __FILE__)
#define _TRACE3() trace3(__LINE__, __FILE__)
#define _TRACE4() trace4(__LINE__, __FILE__)
#define _TRACE5() trace5(__LINE__, __FILE__)


#define CRITICAL_ERROR(__log__, __message__) (__log__)._CRITICAL_ERROR() << __message__ << std::flush << ele
#define ERROR(__log__, __message__) (__log__)._ERROR() << __message__ << std::flush << ele
#define WARNING(__log__, __message__) (__log__)._WARNING() << __message__ << std::flush << ele
#define INFO(__log__, __message__) (__log__)._INFO() << __message__ << std::flush << ele
//#define DEBUG(__log__, __message__) (__log__)._DEBUG() << __message__ << flush << ele
//#define TRACE(__log__, __message__) (__log__)._TRACE() << __message__ << flush << ele
#define LOG_VALUE(__log__, __value_name__, __value__) INFO(__log__, "[" << __value_name__ << "=" << __value__ << "]")
#define LOG_VAL(__log__, __var__)  INFO(__log__, "[" # __var__ "=" << __var__ << "]")


#ifdef LOG_DEBUG
#define DEBUG(__log__, __message__) (__log__)._DEBUG() << __message__ << std::flush << ele
#define TRACE(__log__, __message__)
#define TRACE1(__log__, __message__)
#define TRACE2(__log__, __message__)
#define TRACE3(__log__, __message__)
#define TRACE4(__log__, __message__)
#define TRACE5(__log__, __message__)
#endif

#ifdef LOG_TRACE
#define DEBUG(__log__, __message__) (__log__)._DEBUG() << __message__ << std::flush << ele
#define TRACE(__log__, __message__) (__log__)._TRACE() << __message__ << std::flush << ele
#define TRACE1(__log__, __message__) (__log__)._TRACE1() << __message__ << std::flush << ele
#define TRACE2(__log__, __message__) (__log__)._TRACE2() << __message__ << std::flush << ele
#define TRACE3(__log__, __message__) (__log__)._TRACE3() << __message__ << std::flush << ele
#define TRACE4(__log__, __message__) (__log__)._TRACE4() << __message__ << std::flush << ele
#define TRACE5(__log__, __message__) (__log__)._TRACE5() << __message__ << std::flush << ele
#endif

#if !defined(LOG_TRACE) && !defined(LOG_DEBUG)
#define DEBUG(__log__, __message__)
#define TRACE(__log__, __message__)
#define TRACE1(__log__, __message__)
#define TRACE2(__log__, __message__)
#define TRACE3(__log__, __message__)
#define TRACE4(__log__, __message__)
#define TRACE5(__log__, __message__)
#endif

#define BOOL2STR(b) (b == true ? "true" : "false")


class Logger;
Logger & head(Logger& os);
Logger & ele(Logger & os);


class Logger : public std::basic_ostream<char, std::char_traits<char> > {
  static int log_message_number;

  static log_level_t default_log_level; // = DEFAULT_LOG_LEVEL;
  static bool log_level_read; // = false;
  static int processor_rank_static;
  static int processor_number_static;
  static bool runs_in_parallel;
  //static streambuf *global_stream_buffer;
  static streambuf *global_file_buffer;
  static std::map<string, log_level_t> *log_levels;

  log_level_t loglevel;
  log_level_t _this_entry_loglevel;
  string logger_name;
  string log_head_format;

  //static int global_timestamp;
  int local_timestamp;
  int processor_rank;
  int processor_number;
  streambuf *stream_buffer;

  int line;
  string file;
  bool have_source_information;

  bool _too_detailed;
  bool _writing_le;
  static map<string, Logger *> *loggers;
  void prepare_next_log_entry(int new_line, string &new_file, log_level_t ll) {
    if(_writing_le) (*this) << endl;
    _writing_le = true;
    have_source_information = true;
    this->line = new_line;
    string::size_type tmp = file.find_last_of("/\\");
    if(tmp != string::npos)
      file = file.substr(tmp + 1);
    this->file = new_file;
    _this_entry_loglevel = ll;
    if(loglevel < ll) _too_detailed = true;
  }

  void prepare_next_log_entry(log_level_t ll) {
    if(_writing_le) (*this) << endl;
    _writing_le = true;
    have_source_information = false;
    _this_entry_loglevel = ll;
    if(loglevel < ll) _too_detailed = true;
  }


  void init();

  string process_filename(const string &fname);
public:
  Logger(const Logger &log);

  Logger(string logger_name, const char *head_format, log_level_t loglevel, int default_log_level, const char* filename);
  Logger(string logger_name, const char *head_format, log_level_t loglevel, int default_log_level, const ostream &output);
  Logger(string logger_name, const char *output_format, log_level_t loglevel, int default_log_level, streambuf *obuf);

  virtual ~Logger();

  void set_processor_rank(int pnum) {
    processor_rank = pnum;
  }
  void set_processor_number(int pnum) {
    processor_number = pnum;
  }

  virtual void set_parameter(char *pname, char *value);

  void set_loglevel(log_level_t level) {
    loglevel = level;
  }
  int get_loglevel() {
    return loglevel;
  }

  void set_streambuf(streambuf *buf);

  string get_current_state() {
    std::stringstream ss;
    ss << "_too_detailed=" << _too_detailed << endl;
    ss << "_writing_le=" << _writing_le;
    return ss.str();
  }


  Logger & critical_error() {
    prepare_next_log_entry(LG_CRITICAL_ERROR);
    return (*this << head);
  }

  Logger & critical_error(int l, string f) {
    prepare_next_log_entry(l, f, LG_CRITICAL_ERROR);
    return (*this << head);
  }

  Logger & error() {
    prepare_next_log_entry(LG_ERROR);
    return (*this << head);
  }

  Logger & error(int l, string f) {
    prepare_next_log_entry(l, f, LG_ERROR);
    return (*this << head);
  }

  Logger & warning() {
    prepare_next_log_entry(LG_WARNING);
    return (*this << head);
  }

  Logger & warning(int l, string f) {
    prepare_next_log_entry(l, f, LG_WARNING);
    return (*this << head);
  }

  Logger & info() {
    prepare_next_log_entry(LG_INFO);
    return (*this << head);
  }

  Logger & info(int l, string f) {
    prepare_next_log_entry(l, f, LG_INFO);
    return (*this << head);
  }


  Logger & debug() {
    prepare_next_log_entry(LG_DEBUG);
    return (*this << head);
  }

  Logger &debug(int l, string f) {
    prepare_next_log_entry(l, f, LG_DEBUG);
    return (*this << head);
  }

  Logger & trace() {
    prepare_next_log_entry(LG_TRACE);
    return (*this << head);
  }
  Logger & trace(int l, string f) {
    prepare_next_log_entry(l, f, LG_TRACE);
    return (*this << head);
  }


  Logger & trace1() {
    prepare_next_log_entry(LG_TRACE1);
    return (*this << head);
  }
  Logger & trace1(int l, string f) {
    prepare_next_log_entry(l, f, LG_TRACE1);
    return (*this << head);
  }



  Logger & trace2() {
    prepare_next_log_entry(LG_TRACE1);
    return (*this << head);
  }
  Logger & trace2(int l, string f) {
    prepare_next_log_entry(l, f, LG_TRACE2);
    return (*this << head);
  }



  Logger & trace3() {
    prepare_next_log_entry(LG_TRACE3);
    return (*this << head);
  }
  Logger & trace3(int l, string f) {
    prepare_next_log_entry(l, f, LG_TRACE3);
    return (*this << head);
  }



  Logger & trace4() {
    prepare_next_log_entry(LG_TRACE4);
    return (*this << head);
  }
  Logger & trace4(int l, string f) {
    prepare_next_log_entry(l, f, LG_TRACE4);
    return (*this << head);
  }


  Logger & trace5() {
    prepare_next_log_entry(LG_TRACE5);
    return (*this << head);
  }
  Logger & trace5(int l, string f) {
    prepare_next_log_entry(l, f, LG_TRACE5);
    return (*this << head);
  }



  void too_detailed(bool b = true) {
    _too_detailed = b;
  }
  bool get_detail_flag() {
    return _too_detailed;
  }
  void writing_le(bool b = true) {
    _writing_le = b;
  }

  virtual void set_log_head_format(const char *fmt) {
    log_head_format = (char *) fmt;
  }
  virtual Logger & loghead();


  Logger &operator << (Logger & (*pf)(Logger &)) {
    if(pf == ele) return (*pf)(*this);

    if(_too_detailed) return *this;
    return (*pf)(*this);
  }

  Logger &operator << (Logger & (*pf)(std::basic_ostream<char, std::char_traits<char> > &)) {
    if(_too_detailed) return *this;
    return (*pf)(*this);
  }

  Logger &operator << (std::basic_ostream<char, std::char_traits<char> > & (*pf)(Logger &)) {
    if(_too_detailed) return *this;
    return (Logger &)(*pf)(*this);
  }

  Logger &operator << (std::basic_ostream<char, std::char_traits<char> > & (*pf)(std::basic_ostream<char, std::char_traits<char> > &)) {
    if(_too_detailed) return *this;
    return (Logger &)(*pf)(*this);
  }

  template<class Type> Logger &operator << (Type t) {
    if(_too_detailed) return *this;
    (std::basic_ostream<char, std::char_traits<char> > &) * this << t;
    return *this;
  }

  const char *this_entry_loglevel() {
    return get_log_level(_this_entry_loglevel);
  }

  /*
     static const char *get_log_level(log_level_t l) {
     switch(l) {
     case CRITICAL_ERROR: return "CRITICAL_ERROR";
     case ERROR: return "ERROR";
     case WARNING: return "WARNING";
     case INFO: return "INFO";
     case DEBUG: return "DEBUG";
     case TRACE: return "TRACE";
     case DISABLED: return "DISABLED";
     } // switch
     return "????";
     } // get_log_level
   */
  static const char *get_log_level(log_level_t l) {
    switch(l) {
    case LG_OFF: return "OFF";
    case LG_CRITICAL_ERROR: return "CRITICAL_ERROR";
    case LG_ERROR: return "ERROR";
    case LG_WARNING: return "WARNING";
    case LG_INFO: return "INFO";
    case LG_DEBUG: return "DEBUG";
    case LG_TRACE: return "TRACE";
    case LG_TRACE1: return "TRACE1";
    case LG_TRACE2: return "TRACE2";
    case LG_TRACE3: return "TRACE3";
    case LG_TRACE4: return "TRACE4";
    case LG_TRACE5: return "TRACE5";
    } // switch
      //abort();
    throw std::runtime_error("unknown log level");
  } // get_log_level


  /*
     static log_level_t get_log_level(const char *l) {
     if(strcmp(l, "CRITICAL_ERROR") == 0) return LG_CRITICAL_ERROR;
     else if(strcmp(l, "ERROR") == 0) return LG_ERROR;
     else if(strcmp(l, "WARNING") == 0) return LG_WARNING;
     else if(strcmp(l, "INFO") == 0) return LG_INFO;
     else if(strcmp(l, "DEBUG") == 0) return LG_DEBUG;
     else if(strcmp(l, "TRACE") == 0) return LG_TRACE;
     else if(strcmp(l, "DISABLED") == 0) return LG_DISABLED;
     if(strcmp(l, "CERR") == 0) return LG_CRITICAL_ERROR;
     else if(strcmp(l, "ERR") == 0) return LG_ERROR;
     else if(strcmp(l, "WRN") == 0) return LG_WARNING;
     else if(strcmp(l, "NFO") == 0) return LG_INFO;
     else if(strcmp(l, "DBG") == 0) return LG_DEBUG;
     else if(strcmp(l, "TRC") == 0) return LG_TRACE;
     else if(strcmp(l, "DSBL") == 0) return LG_DISABLED;

     else return DEFAULT_LOG_LEVEL;
     }
   */

  static log_level_t get_log_level(const char *l) {
    if(strcmp(l, "OFF") == 0) return LG_OFF;
    else if(strcmp(l, "LGOFF") == 0) return LG_OFF;
    else if(strcmp(l, "NOLOG") == 0) return LG_OFF;
    else if(strcmp(l, "NL") == 0) return LG_OFF;
    else if(strcmp(l, "CRITICAL_ERROR") == 0) return LG_CRITICAL_ERROR;
    else if(strcmp(l, "CERR") == 0) return LG_CRITICAL_ERROR;
    else if(strcmp(l, "CER") == 0) return LG_CRITICAL_ERROR;
    else if(strcmp(l, "CE") == 0) return LG_CRITICAL_ERROR;
    else if(strcmp(l, "ERROR") == 0) return LG_ERROR;
    else if(strcmp(l, "ERR") == 0) return LG_ERROR;
    else if(strcmp(l, "E") == 0) return LG_ERROR;
    else if(strcmp(l, "WARNING") == 0) return LG_WARNING;
    else if(strcmp(l, "W") == 0) return LG_WARNING;
    else if(strcmp(l, "WRNG") == 0) return LG_WARNING;
    else if(strcmp(l, "INFO") == 0) return LG_INFO;
    else if(strcmp(l, "INF") == 0) return LG_INFO;
    else if(strcmp(l, "I") == 0) return LG_INFO;
    else if(strcmp(l, "NFO") == 0) return LG_INFO;
    else if(strcmp(l, "DEBUG") == 0) return LG_DEBUG;
    else if(strcmp(l, "D") == 0) return LG_DEBUG;
    else if(strcmp(l, "DBG") == 0) return LG_DEBUG;
    else if(strcmp(l, "TRACE") == 0) return LG_TRACE;
    else if(strcmp(l, "T") == 0) return LG_TRACE;
    else if(strcmp(l, "TRC") == 0) return LG_TRACE;
    else if(strcmp(l, "TRACE1") == 0) return LG_TRACE1;
    else if(strcmp(l, "T1") == 0) return LG_TRACE1;
    else if(strcmp(l, "TRC1") == 0) return LG_TRACE1;
    else if(strcmp(l, "TRACE2") == 0) return LG_TRACE2;
    else if(strcmp(l, "T2") == 0) return LG_TRACE2;
    else if(strcmp(l, "TRC2") == 0) return LG_TRACE2;
    else if(strcmp(l, "TRACE3") == 0) return LG_TRACE3;
    else if(strcmp(l, "T3") == 0) return LG_TRACE3;
    else if(strcmp(l, "TRC3") == 0) return LG_TRACE3;
    else if(strcmp(l, "TRACE4") == 0) return LG_TRACE4;
    else if(strcmp(l, "T4") == 0) return LG_TRACE4;
    else if(strcmp(l, "TRC4") == 0) return LG_TRACE4;
    else if(strcmp(l, "TRACE5") == 0) return LG_TRACE5;
    else if(strcmp(l, "T5") == 0) return LG_TRACE5;
    else if(strcmp(l, "TRC5") == 0) return LG_TRACE5;
    else return DEFAULT_LOG_LEVEL;
  }



  static Logger *get_logger(string component);
  static Logger *get_logger(string component, string format);
  static Logger *get_logger(string component, string format, string filename);
  static Logger *get_logger(string component, string format, streambuf *obuf);
  static void free_loggers();
  static void disable_logger(string component);

  static log_level_t get_loglevel_from_env(string logger);

  static void disable_loggers_from_env();
  //static void set_parallel_info(int processor_rank, int processor_number);
  static void set_parallel_info(int processor_rank, int processor_number, bool set_to_existing_loggers = true);
  static void set_parallel_info(int processor_rank, int processor_number, std::string filename_prefix);
  static void set_global_streambuf(streambuf *buf, bool set_to_loggers = false);
};

#endif
