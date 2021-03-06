// Function definitions for MemorySolutionHistory

// Local includes
#include "libmesh/memory_solution_history.h"

#include <cmath>

namespace libMesh
{
// The Destructor
MemorySolutionHistory::~MemorySolutionHistory ()
{
  stored_solutions_iterator stored_sols_it = stored_solutions.begin();
  const stored_solutions_iterator stored_sols_end = stored_solutions.end();

  for (; stored_sols_it != stored_sols_end; ++stored_sols_it)
    {
      // The saved vectors at this timestep
      std::map<std::string, NumericVector<Number> *> saved_vectors = stored_sols_it->second;

      std::map<std::string, NumericVector<Number> *>::iterator vec = saved_vectors.begin();
      std::map<std::string, NumericVector<Number> *>::iterator vec_end = saved_vectors.end();

      // Loop over all the saved vectors
      for (; vec != vec_end; ++vec)
        {
          // Delete this saved vector
          delete vec->second;
        }
    }
}

// This function finds, if it can, the entry where we're supposed to
// be storing data
void MemorySolutionHistory::find_stored_entry()
{
  if (stored_solutions.begin() == stored_solutions.end())
    return;

  libmesh_assert (stored_sols != stored_solutions.end());

  if (std::abs(stored_sols->first - _system.time) < TOLERANCE)
    return;

  // If we're not at the front, check the previous entry
  if (stored_sols != stored_solutions.begin())
    {
      stored_solutions_iterator test_it = stored_sols;
      if (std::abs((--test_it)->first - _system.time) < TOLERANCE)
        {
          --stored_sols;
          return;
        }
    }

  // If we're not at the end, check the subsequent entry
  stored_solutions_iterator test_it = stored_sols;
  if ((++test_it) != stored_solutions.end())
    {
      if (std::abs(test_it->first - _system.time) < TOLERANCE)
        {
          ++stored_sols;
          return;
        }
    }
}

// This functions saves all the 'projection-worthy' system vectors for
// future use
void MemorySolutionHistory::store()
{
  this->find_stored_entry();

  // In an empty history we create the first entry
  if (stored_solutions.begin() == stored_solutions.end())
    {
      stored_solutions.push_back
        (std::make_pair(_system.time,
                        std::map<std::string, NumericVector<Number> *>()));
      stored_sols = stored_solutions.begin();
    }

  // If we're past the end we can create a new entry
  if (_system.time - stored_sols->first > TOLERANCE )
    {
#ifndef NDEBUG
      ++stored_sols;
      libmesh_assert (stored_sols == stored_solutions.end());
#endif
      stored_solutions.push_back
        (std::make_pair(_system.time,
                        std::map<std::string, NumericVector<Number> *>()));
      stored_sols = stored_solutions.end();
      --stored_sols;
    }

  // If we're before the beginning we can create a new entry
  else if (stored_sols->first - _system.time > TOLERANCE)
    {
      libmesh_assert (stored_sols == stored_solutions.begin());
      stored_solutions.push_front
        (std::make_pair(_system.time,
                        std::map<std::string, NumericVector<Number> *>()));
      stored_sols = stored_solutions.begin();
    }

  // We don't support inserting entries elsewhere
  libmesh_assert(std::abs(stored_sols->first - _system.time) < TOLERANCE);

  // Map of stored vectors for this solution step
  std::map<std::string, NumericVector<Number> *> & saved_vectors = stored_sols->second;

  // Loop over all the system vectors
  for (System::vectors_iterator vec = _system.vectors_begin(); vec != _system.vectors_end(); ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      // If we haven't seen this vector before or if we have and
      // want to overwrite it
      if ((overwrite_previously_stored ||
           !saved_vectors.count(vec_name)) &&
          // and if we think it's worth preserving
          _system.vector_preservation(vec_name))
        {
          // Then we save it.
          saved_vectors[vec_name] = vec->second->clone().release();
        }
    }

  // Of course, we will usually save the actual solution
  std::string _solution("_solution");
  if ((overwrite_previously_stored ||
       !saved_vectors.count(_solution)) &&
      // and if we think it's worth preserving
      _system.project_solution_on_reinit())
    saved_vectors[_solution] = _system.solution->clone().release();
}

void MemorySolutionHistory::retrieve()
{
  this->find_stored_entry();

  // Get the time at which we are recovering the solution vectors
  Real recovery_time = stored_sols->first;

  // Print out what time we are recovering vectors at
  //    libMesh::out << "Recovering solution vectors at time: " <<
  //                 recovery_time << std::endl;

  // Do we not have a solution for this time?  Then
  // there's nothing to do.
  if (stored_sols == stored_solutions.end() ||
      std::abs(recovery_time - _system.time) > TOLERANCE)
    {
      //libMesh::out << "No more solutions to recover ! We are at time t = " <<
      //                     _system.time << std::endl;
      return;
    }

  // Get the saved vectors at this timestep
  std::map<std::string, NumericVector<Number> *> & saved_vectors = stored_sols->second;

  std::map<std::string, NumericVector<Number> *>::iterator vec = saved_vectors.begin();
  std::map<std::string, NumericVector<Number> *>::iterator vec_end = saved_vectors.end();

  // Loop over all the saved vectors
  for (; vec != vec_end; ++vec)
    {
      // The name of this vector
      const std::string & vec_name = vec->first;

      // Get the vec_name entry in the saved vectors map and set the
      // current system vec[vec_name] entry to it
      if (vec_name != "_solution")
        _system.get_vector(vec_name) = *(vec->second);
    }

  // Of course, we will *always* have to get the actual solution
  std::string _solution("_solution");
  *(_system.solution) = *(saved_vectors[_solution]);
}

}
// End namespace libMesh
