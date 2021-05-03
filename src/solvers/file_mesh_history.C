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

// C++ include files that we need
#include <iostream>
// Local includes
#include "libmesh/file_mesh_history.h"

#include "libmesh/diff_system.h"
// For the MeshRefinement Object
#include "libmesh/mesh_base.h"
#include "libmesh/mesh.h"
#include "libmesh/inter_mesh_projection.h"
#include "libmesh/equation_systems.h"
#include "libmesh/enum_norm_type.h"

#include <cmath>
#include <iterator>

namespace libMesh
{

  class EquationSystems;
  class MeshBase;

/**
   * Constructor, reference to system to be passed by user, set the
   * stored_meshes_it iterator to some initial value
   */
  FileMeshHistory::FileMeshHistory(System & system_)
  : stored_meshes_it(stored_meshes.end()),
  _system(system_), localTimestamp(0),
  timeTotimestamp()
  {
    libmesh_experimental();
  }


FileMeshHistory::~FileMeshHistory ()
{
}

// This function finds, if it can, the entry where we're supposed to
// be storing data, leaves stored_meshes_it unchanged if it cant find an entry
// with the key corresponding to time.
void FileMeshHistory::find_stored_entry(Real time, bool storing)
{
  if (stored_meshes.begin() == stored_meshes.end())
    return;

  // We will use the map::lower_bound operation to find the key which
  // is the least upper bound among all existing keys for time.
  // (key before map::lower_bound) < time < map::lower_bound, one of these
  // should be within TOLERANCE of time (unless we are creating a new map entry)
  // If the lower bound iterator points to:
  // begin -> we are looking for the mesh at the initial time
  // end -> we are creating a new entry
  // anything else, we are looking for an existing entry
  stored_meshes_iterator lower_bound_it = stored_meshes.lower_bound(time);

  // For the key right before the lower bound
  stored_meshes_iterator lower_bound_it_decremented;

  // If we are at end, we could be creating a new entry (depends on the storing bool), return
  // Otherwise, get a decremented iterator for the sandwich test
  if(lower_bound_it == stored_meshes.end())
  {
    // If we are storing and lower_bound_it points to stored_meshes.end(), we assume
    // that this is a brand new entry in the map. We leave stored_meshes_it unchanged.
    if(storing)
    {
      return;
    }
    else
    {
      // We are trying to retrieve and none of the keys was an upper bound.
      // We could have a situation in which the time is greatest key + FPE.
      // So we can check the key before the end and see if it matches time, else we have an error.
      lower_bound_it = std::prev(lower_bound_it);
    }
  }
  else if(lower_bound_it == stored_meshes.begin()) // At the beginning, so we cant go back any further
  {
    stored_meshes_it = stored_meshes.begin();
    return;
  }
  else // A decremented iterator, to perform the sandwich test for the key closest to time
  {
    lower_bound_it_decremented = std::prev(lower_bound_it);
  }

  // Set the stored sols iterator as per the key which is within TOLERANCE of time
  if(std::abs(lower_bound_it->first - time) < TOLERANCE)
  {
    stored_meshes_it = lower_bound_it;
  }
  else if(std::abs(lower_bound_it_decremented->first - time) < TOLERANCE)
  {
    stored_meshes_it = lower_bound_it_decremented;
  }
  else // Neither of the two candidate keys matched our time
  {
    if(storing) // If we are storing, this is fine, we need to create a new entry, so just return
    {
      return;
    }
    else // If we are not storing, then we expected to find something but didnt, so we have a problem
    {
      libmesh_error_msg("Failed to set stored mesh iterator to a valid value.");
    }
  }


}

// This functions writes the mesh at the current system time to disk
void FileMeshHistory::store(bool is_adjoint_solve, Real time)
{
  // This will map the stored_meshes_it iterator to the current time
  this->find_stored_entry(time, true);

  // In an empty history we create the first entry
  if (stored_meshes.begin() == stored_meshes.end())
    {
      stored_meshes[time] = std::string();
      stored_meshes_it = stored_meshes.begin();
    }

  // If we're past the end we can create a new entry
  if (time - stored_meshes_it->first > TOLERANCE )
    {
#ifndef NDEBUG
      ++stored_meshes_it;
      libmesh_assert (stored_meshes_it == stored_meshes.end());
#endif
      stored_meshes[time] = std::string();
      stored_meshes_it = stored_meshes.end();
      --stored_meshes_it;
    }

  // If we're before the beginning we can create a new entry
  else if (stored_meshes_it->first - time > TOLERANCE)
    {
      libmesh_assert (stored_meshes_it == stored_meshes.begin());
      stored_meshes[time] = std::string();
      stored_meshes_it = stored_meshes.begin();
    }

  // We don't support inserting entries elsewhere
  libmesh_assert(std::abs(stored_meshes_it->first - time) < TOLERANCE);

  // The name of the file to in which we store the mesh from the current timestep
  std::string & mesh_filename = stored_meshes_it->second;

  // This iterator will be used to check if we have already assigned a timestamp for this time key
  // If we have, we are in the adjoint loop, if we have not, we are in the primal loop
  timeTotimestamp_iterator = timeTotimestamp.find(time);

  // Associate the localTimestamp to the current time, if we are in the primal solve loop
  // Then increment the localTimestamp
  if(!is_adjoint_solve)
  {
    // Point mesh_filename to the filename generated by the libMesh I/O object
    mesh_filename = "mesh.out.xda.";
    mesh_filename += std::to_string(localTimestamp);

    // Write the current mesh out to file
    _system.get_mesh().write(mesh_filename);

    timeTotimestamp.insert( std::pair<Real, unsigned int>(time, localTimestamp) );

    ++localTimestamp;
  }
  else // We are in the adjoint time stepping loop
  {
    --localTimestamp;

    // For the adjoint mesh, we reuse the timestamps generated earlier during the primal time march
    mesh_filename = "mesh.out.xda.";
    mesh_filename += std::to_string(localTimestamp);

    // Write the current mesh out to file
    _system.get_mesh().write(mesh_filename);

  }

}

void FileMeshHistory::retrieve(bool is_adjoint_solve, Real time)
{
  this->find_stored_entry(time, false);

  // To set the deltat while using adaptive timestepping, we will utilize
  // consecutive time entries in the stored mesh iterator
  Real _current_time = stored_meshes_it->first;

  // If we are solving the adjoint, we are moving backwards, so decrement time
  // else we are moving forwards, so increment time
  if(is_adjoint_solve)
  {
    stored_meshes_iterator stored_meshes_it_decrement_time = stored_meshes_it;

    // Recovering deltats needs two different entries from the the
    // stored mesh map
    if(stored_meshes_it_decrement_time != stored_meshes.begin())
    {
      stored_meshes_it_decrement_time--;

      Real _decremented_time = stored_meshes_it_decrement_time->first;

      try
      {
        dynamic_cast<DifferentiableSystem &>(_system).deltat = _current_time - _decremented_time;
      }
      catch(const std::bad_cast& e)
      {
        // For a non-diff system, only fixed time step sizes are supported as of now.
      }
    }
  }
  else
  {
    stored_meshes_iterator stored_meshes_it_increment_time = stored_meshes_it;

    // Recovering deltats needs two different entries from the the
    // stored mesh map
    if(stored_meshes_it_increment_time != std::prev(stored_meshes.end()) )
    {
      stored_meshes_it_increment_time++;

      Real _incremented_time = stored_meshes_it_increment_time->first;

      try
      {
        dynamic_cast<DifferentiableSystem &>(_system).deltat = _incremented_time - _current_time;
      }
      catch(const std::bad_cast& e)
      {
        // For a non-diff system, only fixed time step sizes are supported as of now.
      }
    }
  }

    // Get the time at which we are recovering the mesh vectors
    Real recovery_time = stored_meshes_it->first;

    // Do we not have a mesh for this time?  Then
    // there's nothing to do.
    if (stored_meshes_it == stored_meshes.end() || std::abs(recovery_time - time) > TOLERANCE)
    {
      //libMesh::out << "No more mesh to recover ! We are at time t = " <<
      //                     _system.time << std::endl;
      return;
    }

    Real T_H1 = _system.calculate_norm(*_system.solution, 0, H1);
    Real T_L2 = _system.calculate_norm(*_system.solution, 0, L2);

    std::cout << "T_mesh(H1) = " << T_H1 << std::endl << "T_mesh(L2) = " << T_L2 << std::endl;

    // Read in the mesh to be retrieved into a temp_mesh
    Mesh temp_mesh(_system.get_mesh().comm(), _system.get_mesh().mesh_dimension());
    temp_mesh.read(stored_meshes_it->second);

    // Create a temp_system to inter mesh project the solution and system vectors
    // of _system.mesh onto temp_mesh
    EquationSystems temp_mesh_equation_system(temp_mesh);
    System & temp_mesh_system = temp_mesh_equation_system.add_system<System>("temp_system");

    // Add all _system variables to temp_mesh_system
    for ( auto i : make_range(_system.n_vars()) )
      temp_mesh_system.add_variable(_system.variable_name(i), _system.variable(i).type().order, _system.variable(i).type().family);

    temp_mesh_equation_system.init();

    // Add all _system vectors to temp_mesh_system
    for ( System::vectors_iterator vec = _system.vectors_begin(), vec_end = _system.vectors_end(); vec != vec_end; ++vec )
     temp_mesh_system.add_vector(vec->first, temp_mesh_system.vector_preservation(vec->first), (vec->second)->type());

    // Create an IMP object to project from _system to temp_mesh_system
    InterMeshProjection inter_mesh_projector(_system, temp_mesh_system);

    inter_mesh_projector.project_system_vectors();

    Real T_H1_temp_mesh = temp_mesh_system.calculate_norm(*temp_mesh_system.solution, 0, H1);
    Real T_L2_temp_mesh = temp_mesh_system.calculate_norm(*temp_mesh_system.solution, 0, L2);

    std::cout << "T_temp_mesh(H1) = " << T_H1_temp_mesh << std::endl << "T_temp_mesh(L2) = " << T_L2_temp_mesh << std::endl;

    // Reassign _system.mesh and reinit _system's ES with this new mesh
    _system.get_mesh().clear();
    _system.get_mesh().assign(std::move(temp_mesh));
    _system.get_equation_systems().reinit_mesh();

    // Get a pointer to the primal solution vector
    NumericVector<Number> & primal_solution = *_system.solution;

    // Replace the vectors in _system with those from destination_system
    *_system.solution = *temp_mesh_system.solution;

    primal_solution = *_system.solution;

    // Match vectors as well
    for (System::vectors_iterator vec = _system.vectors_begin(), vec_end = _system.vectors_end(); vec != vec_end; ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      _system.get_vector(vec_name) = temp_mesh_system.get_vector(vec_name);
    }

    // We need to call update to put system in a consistent state
    // with the mesh that was read in
    _system.update();

    Real T_H1_assign_mesh = _system.calculate_norm(*_system.solution, 0, H1);
    Real T_L2_assign_mesh = _system.calculate_norm(*_system.solution, 0, L2);

    std::cout << "T_assign_mesh(H1) = " << T_H1_assign_mesh << std::endl << "T_assign_mesh(L2) = " << T_L2_assign_mesh << std::endl;

}

void FileMeshHistory::erase(Real time)
{
  // We cant erase the stored_meshes_it iterator which is used in other places
  // So save its current value for the future
  stored_meshes_iterator stored_meshes_it_last = stored_meshes_it;

  // This will map the stored_meshes_it iterator to the current time
  this->find_stored_entry(time, false);

  // map::erase behaviour is undefined if the iterator is pointing
  // to a non-existent element.
  libmesh_assert(stored_meshes_it != stored_meshes.end());

  // We want to keep using the stored_meshes_it iterator, so we have to create
  // a new one to erase the concerned entry
  stored_meshes_iterator stored_meshes_it_copy = stored_meshes_it;

  // If we're asking to erase the entry at stored_meshes_it, then move stored_meshes_it somewhere safer first
  if(stored_meshes_it == stored_meshes_it_last)
    stored_meshes_it--;

  stored_meshes.erase(stored_meshes_it_copy);

}

}
// End namespace libMesh
