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



#ifndef LIBMESH_PARAMETER_ACCESSOR_H
#define LIBMESH_PARAMETER_ACCESSOR_H


// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"

namespace libMesh
{

// Forward declaration
template <typename T>
class ParameterProxy;


/**
 * Accessor object allowing reading and modification of the
 * independent variables in a parameter sensitivity calculation.
 *
 * This is an abstract base class.  Derived objects may simply modify
 * the parameter value at some address in memory, or may call
 * arbitrary setter/getter functions.
 */
template <typename T=Number>
class ParameterAccessor
{
public:
  /**
   * Setter: change the value of the parameter we access.
   */
  virtual void set (const T & new_value) = 0;

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const T& get () = 0;

  /**
   * Reseater: change the location of the parameter we access.
   * This is included for backward compatibility, but will be
   * deprecated in some classes and not implemented in others.
   */
  virtual void operator= (T * new_ptr)
    { libmesh_error(); }

  /**
   * Proxy: for backward compatibility, we allow codes to treat a
   * ParameterAccessor as if it were a simple pointer-to-value.  We
   * can't safely allow "Number *n = parameter_vector[p]" to compile,
   * but we can allow "*parameter_vector[p] += deltap" to work.
   */

  ParameterProxy<T> operator* ()
    { return ParameterProxy<T>(*this); }
};

template <typename T=Number>
class ParameterProxy
{
public:
  /**
   * Constructor: which parameter are we a proxy for?
   */
  ParameterProxy(ParameterAccessor<T>& accessor) : _accessor(accessor) {}

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy& operator = (const T & new_value)
    { _accessor.set(new_value); }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy& operator += (const T & value_increment)
    { _accessor.set(_accessor.get() + value_increment); }

  /**
   * Setter: change the value of the parameter we access.
   */
  ParameterProxy& operator -= (const T & value_increment)
    { _accessor.set(_accessor.get() - value_increment); }

  /**
   * Getter: get the value of the parameter we access.
   */
  operator T () { return _accessor.get(); }
};


} // namespace libMesh

#endif // LIBMESH_PARAMETER_ACCESSOR_H
