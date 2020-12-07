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

// Local includes
#include "libmesh/inter_mesh_projection.h"

// C++ includes

namespace libMesh
{
    InterMeshProjection::InterMeshProjection(System & _from_system, System & _to_system) :
    from_system(_from_system),
    to_system(_to_system)
    {
        libmesh_experimental();
    }

    void InterMeshProjection::project_system_vectors()
    {
        // Number of vectors to be projected
        unsigned int n_vecs = from_system.n_vectors();
        libmesh_assert_equal_to (to_system.n_vectors(), n_vecs);

        // Number of variable components in each vector
        unsigned int n_vars = from_system.n_vars();
        libmesh_assert_equal_to (to_system.n_vars(), n_vars);

        // Update system vectors with contributions from all processors
        from_system.update();

        // Loop over all the vectors in the system
        for (unsigned int i = 0; i != n_vecs; ++i)
        {
            // A map of mesh functions corresponding to each variable in the from_system
            std::map<std::string, MeshFunction *> mesh_functions;

            // Construct local version of the current system vector
            std::unique_ptr<NumericVector<Number>> current_vector_proxy = NumericVector<Number>::build(from_system.comm());

            std::vector<Number> current_vector;
            current_vector_proxy->init(from_system.get_vector(i).size(), true, SERIAL);
            (*current_vector_proxy) = current_vector;

            // Loop over each variable component of the current system vector
            // for (unsigned int j = 0; j != n_vars; ++j)
            // {
            //     // Construct a MeshFunction for the current component
            //     MeshFunction * mesh_func =
            //     new MeshFunction(from_system.get_equation_systems(), *current_vector_proxy,
            //                  from_system.get_dof_map(), j);

            //     mesh_func->init();
            //     mesh_functions[from_system.variable_name(j)] = mesh_func;
            // }
            // End loop over variable components

            std::vector<unsigned int> variables_vector;

            for (unsigned int j = 0; j != n_vars; ++j)
            {
                variables_vector.push_back(j);
            }

            // Construct a MeshFunction for the current component
            MeshFunction * mesh_func =
            new MeshFunction(from_system.get_equation_systems(), *current_vector_proxy,
                             from_system.get_dof_map(), variables_vector);

            mesh_func->init();

            // Project the current system vector, you need to check if this vector is an adjoint to pass
            // the right options to project_vector
            GradientFunction * gptr = new GradientFunction(mesh_func);
            // FIXME: To properly apply adjoint Dirichlet constraints, we need to pass a fourth argument when we are projecting adjoint vectors.
            to_system.project_vector(to_system.get_vector(i), dynamic_cast<FunctionBase<Number> *>(mesh_func), dynamic_cast<FunctionBase<Gradient> *>(gptr));

            // Delete all the MeshFunctions associated with the current vector
            for (unsigned int j = 0; j != n_vars; ++j)
                delete mesh_functions[from_system.variable_name(j)];
        }
        // End loop over the vectors in the system

    }
    // End InterMeshProjection::project_system_vectors
}
// End namespace libMesh