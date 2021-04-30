// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MEMORY_SOLUTION_HISTORY_H
#define LIBMESH_MEMORY_SOLUTION_HISTORY_H

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/solution_history.h"
#include "libmesh/mesh_history.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ includes
#include <list>

namespace libMesh
{

/**
 * Subclass of Solution History that stores the solutions
 * and other important vectors in memory.
 *
 * \author Vikram Garg
 * \date 2012
 * \brief Stores past solutions in memory.
 */
class MemorySolutionHistory : public SolutionHistory
{
public:

  /**
   * Constructor, reference to system to be passed by user, set the
   * stored_sols iterator to some initial value
   */
  MemorySolutionHistory(System & system_);

  /**
   * Destructor
   */
  ~MemorySolutionHistory();

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
   * Typedef for Stored Solutions iterator, a list of pairs of the current
   * system time, map of strings and saved vectors
   */
  typedef std::map<std::string, std::unique_ptr<NumericVector<Number>>> map_type;
  typedef std::map<Real, map_type> map_map_type;
  typedef map_map_type::iterator stored_solutions_iterator;

  /**
   * Definition of the clone function needed for the setter function
   */
  virtual std::unique_ptr<SolutionHistory > clone() const override
  {
    return libmesh_make_unique<MemorySolutionHistory>(_system);
  }

  /**
   * Local definition of the store_mesh_history() function, which will
   * set the base mesh_history object to type FileMeshHistory and start mesh I/O.
   * This function takes two unsigned ints as argument, to indicate the number of
   * h and p refinements that should be performed on the read in mesh.
   */
  virtual void activate_mesh_history() override;

private:

  // This list of pairs will hold the current time and stored vectors
  // from each timestep
  map_map_type stored_solutions;

  // The stored solutions iterator
  stored_solutions_iterator stored_sols;

  // A helper function to locate entries at a given time
  // Behaviour depends on whether we are calling this function
  // while storing or retrieving/erasing entries.
  // While storing, if no entry in our map matches our time key,
  // we will create a new entry in the map. If we are not storing,
  // not matching a given time key implies an error.
  void find_stored_entry(Real time, bool storing = false);

  // A system reference
  System & _system ;

  /**
   * A vector of pointers to vectors holding the adjoint solution at the last time step
   */
  std::vector< std::unique_ptr<NumericVector<Number>> > dual_solution_copies;

};

} // end namespace libMesh

#endif // LIBMESH_MEMORY_SOLUTION_HISTORY_H
