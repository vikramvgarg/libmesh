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

#ifndef LIBMESH_NO_MESH_HISTORY_H
#define LIBMESH_NO_MESH_HISTORY_H

// Local includes
#include "libmesh/mesh_history.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

/**
 * 'Save nothing' subclass of Mesh History, this is the default.
 *
 * \author Vikram Garg
 * \date 2020
 * \brief For storing and retrieving timestep data.
 */
class NoMeshHistory : public MeshHistory
{
public:

  /**
   * Constructor
   */
  NoMeshHistory() : MeshHistory() {}

  /**
   * Destructor
   */
  virtual ~NoMeshHistory() {}

  /**
   * find_stored_entry will do nothing in this class's implementation.
   */
  virtual void find_stored_entry(Real time, bool storing) override;

  /**
   * Virtual function store which we will be overriding
   */
  virtual void store(bool is_adjoint_solve, Real time) override;

  /**
   * Virtual function retrieve which we will be overriding
   */
  virtual void retrieve(bool is_adjoint_solve, Real time) override;

  /**
   * Virtual function erase which we will be overriding
   */
  virtual void erase(Real time) override;


  /**
   * Definition of the clone function needed for the setter function
   */
  virtual std::unique_ptr<MeshHistory > clone() const
  {
    return libmesh_make_unique<NoMeshHistory>();
  }
};

} // end namespace libMesh

#endif // LIBMESH_NO_MESH_HISTORY_H
