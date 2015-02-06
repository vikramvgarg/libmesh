// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef LIBMESH_PARAMETER_VECTOR_H
#define LIBMESH_PARAMETER_VECTOR_H


// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_pointer.h"

// C++ Includes   -----------------------------------
#include <vector>

namespace libMesh
{


/**
 * Data structure for specifying which Parameters should be
 * independent variables in a parameter sensitivity calculation.
 */
class ParameterVector
{
public:
  /**
   * Default constructor: "no parameters"
   */
  ParameterVector() : _is_shallow_copy(false) {}

  /**
   * Constructor-from-vector-of-Number*: each points to a parameter
   */
  explicit
  ParameterVector(const std::vector<Number *> &params);

  /**
   * Destructor - deletes ParameterAccessor objects
   */
  ~ParameterVector();

  /**
   * Deep copy constructor: the \p target will now own new copies of
   * all the parameter values I'm pointing to
   */
  void deep_copy(ParameterVector &target) const;

  /**
   * Shallow copy constructor: the \p target will now point to all the
   * parameter values I'm pointing to
   */
  void shallow_copy(ParameterVector &target) const;

  /**
   * Value copy method: the \p target, which should already have as
   * many parameters as I do, will now have those parameters set to my
   * values.
   */
  void value_copy(ParameterVector &target) const;

  /**
   * Resets to "no parameters"
   */
  void clear();

  /**
   * Returns the number of parameters to be used
   */
  std::size_t size() const { return _params.size(); }

  /**
   * Sets the number of parameters to be used.  This method is for
   * resizing a ParameterVector that acts as a proxy to other
   * parameter values
   */
  void resize(unsigned int s) { _params.resize(s); }

  /**
   * Sets the number of parameters to be used.  This method is for
   * resizing a ParameterVector that owns its own parameter values
   */
  void deep_resize(unsigned int s);

  /**
   * Returns a smart-pointer to a parameter value
   */
  const ParameterAccessor<Number>& operator[](unsigned int i) const;

  /**
   * Returns a reference to a smart-pointer to a parameter value,
   * suitable for repointing it to a different address.
   * This method is deprecated and may not work with more
   * sophisticated ParameterAccessor subclasses.
   */
  ParameterAccessor<Number>& operator[](unsigned int i);

  /**
   * Multiplication operator; acts individually on each parameter.
   */
  ParameterVector& operator *= (const Number a);

  /**
   * Addition operator.  The parameter vector to be added in must
   * have the same number of values.
   */
  ParameterVector& operator += (const ParameterVector& a);

private:
  /**
   * Pointers to parameters which may exist elsewhere
   */
  std::vector<ParameterAccessor<Number> *> _params;

  /**
   * Parameters which I own; e.g. as the result of a deep copy
   */
  std::vector<Number> _my_data;

  /**
   * Am I a shallow copy?  If so then I shouldn't be deleting my
   * ParameterAccessors.
   */
  bool _is_shallow_copy;
};



// ------------------------------------------------------------
// ParameterVector inline methods

inline
ParameterVector::ParameterVector(const std::vector<Number *> &params)
  : _is_shallow_copy(false)
{
  _params.reserve(params.size());

  for (unsigned int i=0; i != params.size(); ++i)
    _params.push_back(new ParameterPointer<Number>(params[i]));
}


inline
ParameterVector::~ParameterVector()
{
  this->clear();
}


inline
void
ParameterVector::clear()
{
  if (_is_shallow_copy)
    for (unsigned int i=0; i != _params.size(); ++i)
      delete _params[i];

  _params.clear();
  _my_data.clear();
}



inline
const ParameterAccessor<Number>& ParameterVector::operator[] (unsigned int i) const
{
  libmesh_assert_greater (_params.size(), i);

  return *_params[i];
}



inline
ParameterAccessor<Number>& ParameterVector::operator[] (unsigned int i)
{
  libmesh_assert_greater (_params.size(), i);

  return *_params[i];
}

} // namespace libMesh

#endif // LIBMESH_PARAMETER_VECTOR_H
