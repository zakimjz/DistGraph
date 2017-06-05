/*
 * This source is part of the single graph mining algorithm.
 *
 * Copyright (C) 2014 Robert Kessl
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __SERIALIZABLE_HPP__
#define __SERIALIZABLE_HPP__

namespace types {

class serializable_buffer {
public:
  virtual size_t get_serialized_size() const = 0;
  virtual size_t get_serialized_size(char *buffer, size_t buffer_size) const = 0;
  virtual size_t serialize(char *buffer, size_t buffer_size) const = 0;
  virtual size_t deserialize(char *buffer, size_t buffer_size) = 0;
};


class serializable_stream {
public:
  virtual size_t serialize(std::ostream &) const = 0;
  virtual size_t deserialize(std::istream &) = 0;
};

class serializable : public serializable_buffer, public serializable_stream {
public:
  virtual ~serializable() {
  }
};


} // namespace types

#endif

