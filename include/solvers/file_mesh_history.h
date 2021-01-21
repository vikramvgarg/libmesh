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



#ifndef LIBMESH_FILE_MESH_HISTORY_H
#define LIBMESH_FILE_MESH_HISTORY_H

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_history.h"
#include "libmesh/distributed_mesh.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/inter_mesh_projection.h"
#include "libmesh/libmesh.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <list>

namespace libMesh
{

/**
 * Subclass of Mesh History that stores the mesh at each timestep onto the disk.
 *
 * \author Vikram Garg
 * \date 2020
 * \brief Stores past meshes onto disk.
 */
class FileMeshHistory : public MeshHistory
{
public:

  /**
   * Constructor, reference to system to be passed by user, set the
   * stored_sols iterator to some initial value
   */
  FileMeshHistory(System & system_);

  /**
   * Destructor
   */
  ~FileMeshHistory();

  /**
   * Virtual function store which we will be overriding to store timesteps
   */
  virtual void store(bool is_adjoint_solve, Real time) override;

  /**
   * Virtual function retrieve which we will be overriding to retrieve timesteps
   */
  virtual void retrieve(bool is_adjoint_solve, Real time) override;

  /**
   * Virtual function retrieve which we will be overriding to erase timesteps
   */
  virtual void erase(Real time) override;

  /**
   * Typedef for Stored Meshes iterator, a list of pairs of the
   * system time and filenames of the stored meshes
   */
  typedef std::map<Real, std::string> map_type;
  typedef map_type::iterator stored_meshes_iterator;

  /**
   * Definition of the clone function needed for the setter function
   */
  virtual std::unique_ptr<MeshHistory > clone() const override
  {
    return libmesh_make_unique<FileMeshHistory>(_system);
  }

private:

  // This list of pairs will hold the timestamp and filename of each stored mesh
  map_type stored_meshes;

  // The stored meshes iterator
  stored_meshes_iterator stored_meshes_it;

  // A helper function to locate entries at a given time
  // Behaviour depends on whether we are calling this function
  // while storing or retrieving/erasing entries.
  // While storing, if no entry in our map matches our time key,
  // we will create a new entry in the map. If we are not storing,
  // not matching a given time key implies an error.
  void find_stored_entry(Real time, bool storing = false);

  // A system reference
  System & _system ;

  // A 'timestamp' that belongs specifically to FMH, this will be used to generate filenames
  unsigned int localTimestamp;

  // To assign filenames a timestamp, we will maintain a datastructure within
  // FileMeshHistory which will map system.time to 'timestamps'
  std::map<Real, unsigned int> timeTotimestamp;

  // An iterator for the timeTotimestamp map, will help member functions distinguish
  // between primal and adjoint time loops
  std::map<Real, unsigned int>::iterator timeTotimestamp_iterator;

  // ES & System reference to hold the projected vectors provided by inter_mesh_projector
  //MeshBase * target_mesh;
  //EquationSystems * inter_mesh_equation_systems;
  //System * inter_mesh_system;

  // An inter mesh projection object to project extant adjoint and other system vectors
  // onto the incoming mesh
  //InterMeshProjection * inter_mesh_projector;
};

} // end namespace libMesh

#endif // LIBMESH_FILE_MESH_HISTORY_H
