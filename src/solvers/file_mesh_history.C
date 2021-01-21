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
#include "libmesh/boundary_info.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system_norm.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/dof_map.h"

#include <cmath>
#include <iterator>

namespace libMesh
{

  class EquationSystems;
  class MeshBase;
  class BoundaryInfo;

/**
   * Constructor, reference to system to be passed by user, set the
   * stored_meshes_it iterator to some initial value
   */
  FileMeshHistory::FileMeshHistory(System & system_)
  : stored_meshes_it(stored_meshes.end()),
  _system(system_), localTimestamp(0),
  timeTotimestamp() //, target_mesh(&system_.get_mesh())
  {
    // Intialize the inter_mesh_equation_system with a mesh of the right dimension
    //target_mesh = new Mesh(system_.get_mesh().comm(), system_.get_mesh().mesh_dimension());
    //inter_mesh_equation_systems = new EquationSystems(*target_mesh);

    // Set up the inter_mesh_system to be supplied to inter_mesh_projector with the
    // the right variables and vectors
    //inter_mesh_system = & inter_mesh_equation_systems->add_system<System> ("InterMeshSystem");

    // Match variables
    //unsigned int n_vars = _system.n_vars();
    //for (unsigned int j = 0; j != n_vars; ++j)
    //{
      //inter_mesh_system->add_variable(_system.variable_name(j), _system.variable(j).type().order, _system.variable(j).type().family);
    //}

    // Match vectors as well
    //for (System::vectors_iterator vec = _system.vectors_begin(), vec_end = _system.vectors_end(); vec != vec_end; ++vec)
    //{
      //inter_mesh_system->add_vector(vec->first, inter_mesh_system->vector_preservation(vec->first), (vec->second)->type());
    //}

    //inter_mesh_equation_systems->init ();

    //inter_mesh_projector = new InterMeshProjection(_system, *inter_mesh_system);

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

    // Testing norms of system vectors
    // Compute the H1 norm of the velocity variables
    Real T_H1 = _system.calculate_norm(*_system.solution, 0, H1);
    Real T_L2 = _system.calculate_norm(*_system.solution, 0, L2);

    std::cout << "T_before(H1) = " << T_H1 << std::endl << "T_before(L2) = " << T_L2 << std::endl;

    // Compute the H1 norm of the velocity variables
	  Real T_adjoint_H1 = _system.calculate_norm(_system.get_adjoint_solution(), 0, H1);
	  Real T_adjoint_L2 = _system.calculate_norm(_system.get_adjoint_solution(), 0, L2);

	  std::cout << "T_adjoint(H1) = " << T_adjoint_H1 << std::endl << "T_adjoint(L2) = " << T_adjoint_L2 << std::endl ;

    // Print solution vector
    libMesh::out << "U = [" << (*_system.solution)
                   << "];" << std::endl;

    // Print adjoint vector
    libMesh::out << "Z = [" << _system.get_adjoint_solution()
                   << "];" << std::endl;

    // Testing that the mesh properties match
    _system.get_mesh().print_info();

    Mesh temp_mesh(_system.get_mesh().comm(), _system.get_mesh().mesh_dimension());
    temp_mesh.read(stored_meshes_it->second);
    temp_mesh.print_info();

    EquationSystems inter_mesh_equation_systems(temp_mesh);
    //EquationSystems inter_mesh_equation_systems(_system.get_mesh());
    System & inter_mesh_system = inter_mesh_equation_systems.add_system<System> ("temp_system");
    // Match variables
    unsigned int n_vars = _system.n_vars();
    for (unsigned int j = 0; j != n_vars; ++j)
    {
      inter_mesh_system.add_variable(_system.variable_name(j), _system.variable(j).type().order, _system.variable(j).type().family);
    }

    inter_mesh_equation_systems.init();

    // Match vectors as well
    for (System::vectors_iterator vec = _system.vectors_begin(), vec_end = _system.vectors_end(); vec != vec_end; ++vec)
    {
      inter_mesh_system.add_vector(vec->first, inter_mesh_system.vector_preservation(vec->first), (vec->second)->type());
    }
    InterMeshProjection inter_mesh_projector(_system, inter_mesh_system);

    inter_mesh_projector.project_system_vectors();

    // Testing norms of system vectors
    // Compute the H1 norm of the velocity variables
    T_H1 = _system.calculate_norm(*_system.solution, 0, H1);
    T_L2 = _system.calculate_norm(*_system.solution, 0, L2);

    std::cout << "T_after_projection(H1) = " << T_H1 << std::endl << "T_after_projection(L2) = " << T_L2 << std::endl;

    T_adjoint_H1 = _system.calculate_norm(_system.get_adjoint_solution(), 0, H1);
	  T_adjoint_L2 = _system.calculate_norm(_system.get_adjoint_solution(), 0, L2);

	  std::cout << "T_adjoint_after_projection(H1) = " << T_adjoint_H1 << std::endl << "T_adjoint_after_projection(L2) = " << T_adjoint_L2 << std::endl ;

    // // At this point, we could replaced the old mesh with the new read in mesh, but then we wont be able to call es.reinit() to project
    // // vectors, since the projections there are meant to be performed between different refinement levels
    // // of the same mesh.
    // // Instead we will use an InterMeshProjection object and a System _projection_system object
    // // to project our vectors from the old mesh to the new one being read in.
    // // Then, we can clear the old mesh, read in the new one and replace
    // // vectors in _system with those from _projection_system and carry on.
    //Mesh temp_mesh(_system.get_mesh().comm(), _system.get_mesh().mesh_dimension());
    //temp_mesh.read(stored_meshes_it->second);
    //inter_mesh_system->get_mesh().clear();
    //inter_mesh_system->get_mesh().assign(temp_mesh);

    // Before calling the IMP methods, we need to match the variables and vectors
    // between _system and inter_mesh_system

    //inter_mesh_equation_systems->reinit();
    // EquationSystems destination_equation_systems (destination_mesh);
    // System &destination_system = destination_equation_systems.add_system<System> ("DestinationHeatSystem");
    // // Add all the variables from _system to destination_system
    // // Add all the vectors from _system to destination_system
    //InterMeshProjection inter_mesh_projector(_system , *inter_mesh_system);
    //inter_mesh_projector->project_system_vectors();

    // Read in the mesh at this time instant and reinit to project solutions on to the new mesh
    _system.get_mesh().clear();
    //_system.get_mesh().read(stored_meshes_it->second);
    _system.get_mesh().assign(temp_mesh);
    _system.get_mesh().print_info();

    // Reinit the DofMap now that the mesh has changed
    std::cout<<"N_systems_1: "<<_system.get_equation_systems().n_systems()<<std::endl;
    //_system.get_dof_map().distribute_dofs(_system.get_mesh());
    _system.get_equation_systems().init();

    // Replace the vectors in _system with those from destination_system
    _system.solution = std::move(inter_mesh_system.solution);

    // Match vectors as well
    for (System::vectors_iterator vec = _system.vectors_begin(), vec_end = _system.vectors_end(); vec != vec_end; ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      _system.get_vector(vec_name) = inter_mesh_system.get_vector(vec_name);
    }

    // // Testing norms of system vectors
    // // Compute the H1 norm of the velocity variables
    // T_H1 = _system.calculate_norm(*_system.solution, 0, H1);
    // T_L2 = _system.calculate_norm(*_system.solution, 0, L2);

    // std::cout << "T_after_assign(H1) = " << T_H1 << std::endl << "T_after_assign(L2) = " << T_L2 << std::endl;

    // T_adjoint_H1 = _system.calculate_norm(_system.get_adjoint_solution(), 0, H1);
	  // T_adjoint_L2 = _system.calculate_norm(_system.get_adjoint_solution(), 0, L2);

	  // std::cout << "T_adjoint_after_assign(H1) = " << T_adjoint_H1 << std::endl << "T_adjoint_after_assign(L2) = " << T_adjoint_L2 << std::endl ;

    // // Print solution vector
    // libMesh::out << "U_after_assign = [" << (*_system.solution)
    //                << "];" << std::endl;

    // // Print adjoint vector
    // libMesh::out << "Z_after_assign = [" << _system.get_adjoint_solution()
    //                << "];" << std::endl;

    std::cout<<"N_systems_2: "<<_system.get_equation_systems().n_systems()<<std::endl;
    _system.get_equation_systems().reinit();

    // We need to call update to put system in a consistent state
    // with the mesh that was read in
    _system.update();

    // Testing norms of system vectors
    // Compute the H1 norm of the velocity variables
    T_H1 = _system.calculate_norm(*_system.solution, 0, H1);
    T_L2 = _system.calculate_norm(*_system.solution, 0, L2);

    std::cout << "T_after_reinit(H1) = " << T_H1 << std::endl << "T_after_assign(L2) = " << T_L2 << std::endl;

    T_adjoint_H1 = _system.calculate_norm(_system.get_adjoint_solution(), 0, H1);
	  T_adjoint_L2 = _system.calculate_norm(_system.get_adjoint_solution(), 0, L2);

	  std::cout << "T_adjoint_after_reinit(H1) = " << T_adjoint_H1 << std::endl << "T_adjoint_after_assign(L2) = " << T_adjoint_L2 << std::endl ;

    // Print solution vector
    libMesh::out << "U_after_reinit = [" << (*_system.solution)
                   << "];" << std::endl;

    // Print adjoint vector
    libMesh::out << "Z_after_reinit = [" << _system.get_adjoint_solution()
                   << "];" << std::endl;

    // Testing
    auto n_conds = _system.get_mesh().get_boundary_info().n_boundary_conds();

    std::cout<<"n_conds: "<<n_conds;

    _system.get_mesh().write("test.xda");

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
