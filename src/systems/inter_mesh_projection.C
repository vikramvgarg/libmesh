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

#include "libmesh/system_norm.h"
#include "libmesh/enum_norm_type.h"


// C++ includes

namespace libMesh
{
    GradientFunction::GradientFunction(MeshFunction * _mesh_function)
    {
     //mesh_function = std::unique_ptr<MeshFunction>(dynamic_cast<MeshFunction *>((_mesh_function->clone()).get()));
     mesh_function = _mesh_function->clone();
     libmesh_experimental();
    }

    void GradientFunction::operator() (const Point & p, const Real, DenseVector<Gradient> & output)
    {
        Real time = 0.0;
        (dynamic_cast<MeshFunction *>(mesh_function.get()))->gradient(p, time, output.get_values());
        return;
    }

    InterMeshProjection::InterMeshProjection(System & _from_system, System & _to_system) :
    from_system(_from_system),
    to_system(_to_system)
    {
        libmesh_experimental();
    }

    void InterMeshProjection::project_system_vectors()
    {
        libmesh_assert_equal_to (to_system.n_vectors(), from_system.n_vectors() );

        // Number of variable components in each vector
        unsigned int n_vars = from_system.n_vars();
        libmesh_assert_equal_to (to_system.n_vars(), n_vars);

        // We are going to use the multi-variable MeshFunction, so we can pass
        // a single vector of variables rather than have a MeshFunction for each variable
	    std::vector<unsigned int> variables_vector;

        for (unsigned int j = 0; j != n_vars; ++j)
        {
          variables_vector.push_back(j);
        }

        // Update system vectors with contributions from all processors
        from_system.update();

        // Any system holds the solution along with the other vectors system.vectors
        // We will first project the solution and then move to the system.vectors

        // Construct local version of the current system vector
	    // This has to be a serial vector
        std::unique_ptr<NumericVector<Number>> solution_vector_serial = NumericVector<Number>::build(from_system.comm());
	    solution_vector_serial->init(from_system.solution->size(), true, SERIAL);

        std::vector<Number> solution_vector;
        from_system.update_global_solution(solution_vector);
        (*solution_vector_serial) = solution_vector;

        // Construct a MeshFunction for the solution
        MeshFunction * mesh_func_solution =
            new MeshFunction(from_system.get_equation_systems(), *solution_vector_serial,
                             from_system.get_dof_map(), variables_vector);

        mesh_func_solution->init();

        // For some element types (say C1) we also need to pass a gradient evaluation MeshFunction
        // To do this evaluate, a new shim class GradientMeshFunction has been added which redirects
        // gptr::operator evaluations inside projection methods into MeshFunction::gradient calls.
        GradientFunction * gptr_solution = new GradientFunction(mesh_func_solution);
        gptr_solution->init();

        from_system.update();

        // // Testing
        // Point p_eval_sol;
        // p_eval_sol(0) = 0.6;
        // //p_eval(1) = 0.5;

        // DenseVector<Number> vector_point_value_sol;
        // mesh_func_solution->operator()(p_eval_sol, 0.0, vector_point_value_sol);

        // std::cout<<"vector_point_value_sol(0.6): "<<vector_point_value_sol(0)<<std::endl;

        // DenseVector<Gradient> vector_gradient_from_mesh_func_sol;
        // mesh_func_solution->gradient(p_eval_sol, 0.0, vector_gradient_from_mesh_func_sol.get_values());
        // std::cout<<"gradient_point_value_from_mesh_func_sol(0.6): "<<vector_gradient_from_mesh_func_sol(0)(0)<<std::endl;

        // DenseVector<libMesh::Gradient> vector_gradient_from_gptr_sol;
        // gptr_solution->operator()(p_eval_sol, 0.0, vector_gradient_from_gptr_sol);
        // std::cout<<"gradient_point_value_sol(0.6): "<<vector_gradient_from_gptr_sol(0)(0)<<std::endl;

        Real T_from_H1 = from_system.calculate_norm(*from_system.solution, 0, H1);
        Real T_from_L2 = from_system.calculate_norm(*from_system.solution, 0, L2);

        // Print the vector
        for(unsigned int i = 0; i < from_system.solution->size(); i++)
            std::cout<<(*from_system.solution)(i)<<std::endl;

        std::cout << "T_from(H1) = " << T_from_H1 << std::endl << "T_from(L2) = " << T_from_L2 << std::endl;

        // End testing

        // FIXME: To properly apply adjoint Dirichlet constraints, we need to pass a fourth argument when we are projecting adjoint vectors.
        to_system.project_vector(*to_system.solution, dynamic_cast<FunctionBase<Number> *>(mesh_func_solution), dynamic_cast<FunctionBase<Gradient> *>(gptr_solution));

        //to_system.get_equation_systems().reinit();

        // Begin testing
        Real T_to_H1 = to_system.calculate_norm(*to_system.solution, 0, H1);
        Real T_to_L2 = to_system.calculate_norm(*to_system.solution, 0, L2);

        // Print the vector
        for(unsigned int i = 0; i < to_system.solution->size(); i++)
            std::cout<<(*to_system.solution)(i)<<std::endl;

        std::cout << "T_to(H1) = " << T_to_H1 << std::endl << "T_to(L2) = " << T_to_L2 << std::endl;

        // Construct local version of the current system vector
        std::unique_ptr<NumericVector<Number>> to_solution_proxy = NumericVector<Number>::build(from_system.comm());

        std::vector<Number> to_solution;
        to_system.update_global_solution(to_solution);
        to_solution_proxy->init(to_system.solution->size(), true, SERIAL);
        (*to_solution_proxy) = to_solution;

        // Construct a MeshFunction for the current component
        MeshFunction * to_mesh_func_sol =
            new MeshFunction(to_system.get_equation_systems(), *to_solution_proxy,
                             to_system.get_dof_map(), variables_vector);

        to_mesh_func_sol->init();

        // DenseVector<Number> vector_point_value_to_mesh_func_sol;
        // to_mesh_func_sol->operator()(p_eval_sol, 0.0, vector_point_value_to_mesh_func_sol);

        // std::cout<<"vector_point_value_sol(0.6): "<<vector_point_value_to_mesh_func_sol(0)<<std::endl;

        // DenseVector<Gradient> vector_gradient_to_mesh_func_sol;
        // to_mesh_func_sol->gradient(p_eval_sol, 0.0, vector_gradient_to_mesh_func_sol.get_values());
        // std::cout<<"gradient_point_value_to_mesh_func_sol(0.6): "<<vector_gradient_to_mesh_func_sol(0)(0)<<std::endl;

        // End testing

        // Delete all the MeshFunctions associated with the solution
        delete mesh_func_solution;
        delete gptr_solution;
        delete to_mesh_func_sol;

        // Now loop over the vectors in system.vectors (includes old_nonlin_sol, rhs, adjoints, adjoint_rhs, sensitivity_rhs)
    	for (System::vectors_iterator vec = from_system.vectors_begin(), vec_end = from_system.vectors_end(); vec != vec_end; ++vec)
	    {
	        // The name of this vector
	        const std::string & vec_name = vec->first;

            std::cout<<"Projecting vector named: "<<vec_name<<std::endl;

            // Print the vector
            for(unsigned int i = 0; i < from_system.get_vector(vec_name).size(); i++)
                std::cout<<(from_system.get_vector(vec_name))(i)<<std::endl;

            // Construct local version of the current system vector
	        // This has to be a serial vector
            std::unique_ptr<NumericVector<Number>> current_vector_proxy = NumericVector<Number>::build(from_system.comm());
	        current_vector_proxy->init(from_system.get_vector(vec_name).size(), true, SERIAL);

            from_system.get_vector(vec_name).localize(*current_vector_proxy);

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

            // Construct a MeshFunction for the current component
            MeshFunction * mesh_func =
            new MeshFunction(from_system.get_equation_systems(), *current_vector_proxy,
                             from_system.get_dof_map(), variables_vector);

            mesh_func->init();

            // Project the current system vector, you need to check if this vector is an adjoint to pass
            // the right options to project_vector
            GradientFunction * gptr = new GradientFunction(mesh_func);
            gptr->init();

            from_system.update();

            // Begin testing
            Point p_eval_vec;
            p_eval_vec(0) = 0.5;
            //p_eval(1) = 0.5;

            DenseVector<Number> vector_point_value;
            mesh_func->operator()(p_eval_vec, 0.0, vector_point_value);

            std::cout<<"vector_point_value(0.5): "<<vector_point_value(0)<<std::endl;

            DenseVector<Gradient> vector_gradient_from_mesh_func;
            mesh_func->gradient(p_eval_vec, 0.0, vector_gradient_from_mesh_func.get_values());
            std::cout<<"gradient_point_value_from_mesh_func(0.5): "<<vector_gradient_from_mesh_func(0)(0)<<std::endl;

            DenseVector<libMesh::Gradient> vector_gradient_from_gptr;
            gptr->operator()(p_eval_vec, 0.0, vector_gradient_from_gptr);
            std::cout<<"gradient_point_value(0.5): "<<vector_gradient_from_gptr(0)(0)<<std::endl;

            Real vec_from_H1 = from_system.calculate_norm(from_system.get_vector(vec_name), 0, H1);
            Real vec_from_L2 = from_system.calculate_norm(from_system.get_vector(vec_name), 0, L2);

            std::cout << "vec_from(H1) = " << vec_from_H1 << std::endl << "vec_from(L2) = " << vec_from_L2 << std::endl;

            // FIXME: To properly apply adjoint Dirichlet constraints, we need to pass a fourth argument when we are projecting adjoint vectors.
            to_system.project_vector(to_system.get_vector(vec_name), dynamic_cast<FunctionBase<Number> *>(mesh_func), dynamic_cast<FunctionBase<Gradient> *>(gptr));

            to_system.update();

            // Print the vector
            for(unsigned int i = 0; i < to_system.get_vector(vec_name).size(); i++)
                std::cout<<(to_system.get_vector(vec_name))(i)<<std::endl;

            Real vec_to_H1 = to_system.calculate_norm(to_system.get_vector(vec_name), 0, H1);
            Real vec_to_L2 = to_system.calculate_norm(to_system.get_vector(vec_name), 0, L2);

            std::cout << "vec_to(H1) = " << vec_to_H1 << std::endl << "vec_to(L2) = " << vec_to_L2 << std::endl<<std::endl;

            // Construct local version of the current system vector
            std::unique_ptr<NumericVector<Number>> to_current_vector_proxy = NumericVector<Number>::build(from_system.comm());

            to_current_vector_proxy->init(to_system.get_vector(vec_name).size(), true, SERIAL);
            to_system.get_vector(vec_name).localize(*to_current_vector_proxy);

            // Construct a MeshFunction for the current component
            MeshFunction * to_mesh_func =
            new MeshFunction(to_system.get_equation_systems(), *to_current_vector_proxy,
                             to_system.get_dof_map(), variables_vector);

            to_mesh_func->init();

            DenseVector<Number> vector_point_value_to_mesh_func;
            to_mesh_func->operator()(p_eval_vec, 0.0, vector_point_value_to_mesh_func);
            std::cout<<"vector_point_value(0.5): "<<vector_point_value_to_mesh_func(0)<<std::endl;

            DenseVector<Gradient> vector_gradient_to_mesh_func;
            to_mesh_func->gradient(p_eval_vec, 0.0, vector_gradient_to_mesh_func.get_values());
            std::cout<<"gradient_point_value_to_mesh_func(0.5): "<<vector_gradient_to_mesh_func(0)(0)<<std::endl;

            // Delete all the MeshFunctions associated with the current vector
            delete mesh_func;
            delete gptr;
            delete to_mesh_func;
        }
        // End loop over the vectors in the system

    }
    // End InterMeshProjection::project_system_vectors
}
// End namespace libMesh