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



#ifndef LIBMESH_INTER_MESH_PROJECTION_H
#define LIBMESH_INTER_MESH_PROJECTION_H

// Local includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/system.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/fem_function_base.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{

    // Forward declarations

    /**
     * This class implements inter mesh projection, i.e. projection of
     * vectors defined on a given mesh (from_mesh associated with from_system)
     * to another mesh (to_mesh), stored in the map projected_vectors_map.
     */

    class InterMeshProjection
    {
        public:

        // Constructor, specifies the _from_system whose vectors will be
        // projected onto the _to_mesh
        InterMeshProjection(System & _from_system, System & _to_mesh);

        // Projects from_system vectors onto the to_mesh
        void project_system_vectors();

        static Number fptr(const Point & p, const Parameters &, const std::string & libmesh_dbg_var(sys_name), const std::string & unknown_name);

        static Gradient gptr(const Point & p, const Parameters &, const std::string & libmesh_dbg_var(sys_name), const std::string & unknown_name);

        private:

        // Local copy of the _from_system
        System & from_system;

        // Local copy of the _to_system
        System & to_system;

    };

    class GradientFunction : public FunctionBase<Gradient>
    {
        public:
        // Constructor
        GradientFunction(MeshFunction * _mesh_function)
        {mesh_function = _mesh_function;}

        // Destructor
        virtual ~GradientFunction () {}

        virtual std::unique_ptr<FunctionBase<Gradient>> clone () const
        {
            return libmesh_make_unique<GradientFunction>(*this);
        }

        virtual Gradient operator() (const Point & ,
                             const Real)
        { libmesh_not_implemented(); }

        virtual void operator() (const Point & p,
                             const Real,
                             DenseVector<Gradient> & output)
        { Real time = 0.0; mesh_function->gradient(p, time, dynamic_cast<std::vector<Gradient> &>(output)); return;}

        private:
        MeshFunction * mesh_function;

    };
}

#endif // LIBMESH_INTER_MESH_PROJECTION_H