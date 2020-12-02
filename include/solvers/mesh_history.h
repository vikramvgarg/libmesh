// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_MESH_HISTORY_H
#define LIBMESH_MESH_HISTORY_H

// Local Includes
#include "libmesh/system.h"

namespace libMesh
{

/**
 * A MeshHistory class that enables the storage and retrieval of
 * the mesh at timesteps.
 *
 * \author Vikram Garg
 * \date 2020
 * \brief For storing and retrieving timestep mesh data.
 */
class MeshHistory
{
public:

  /**
   * Constructor
   */
  MeshHistory(){}

  /**
   * Destructor
   */
  virtual ~MeshHistory () {}

  /**
   * Function to find a stored entry, pure virtual
   */
  virtual void find_stored_entry(Real /*time*/, bool /*storing*/) = 0;

  /**
   * Function to store a solution
   */
  virtual void store(bool /*is_adjoint_solve*/, Real /*time*/) = 0;

  /**
   * Function to retrieve a solution
   */
  virtual void retrieve(bool /*is_adjoint_solve*/, Real /*time*/) = 0;

  /**
   * Function to erase solution at a given time
   */
  virtual void erase(Real /*time*/) = 0;

  /**
   * Cloning function for a std::unique_ptr.
   */
  virtual std::unique_ptr<MeshHistory > clone() const = 0;

};

} // end namespace libMesh

#endif // LIBMESH_MESH_HISTORY_H
