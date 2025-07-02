# Copyright (C) 2022 TurboRVB group
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
This script generates fortran wrapper for cudalib C API

Ot(t)o Kohu\'ak
"""

import csnake as cs
from includer import add_this

cw = cs.CodeWriter()

from basics import ( add_complex,
                     default_includes,
                     default_definitions,
                     add_notice,
                     add_license,
                   )

add_license(cw)
add_notice(cw)

default_includes(cw)
default_definitions(cw)

add_complex(cw)

b_cublas = False
b_cublas_v2 = True
b_cusolver = True

if b_cublas and b_cublas_v2:
    raise RuntimeError("b_cublas and b_cublas_v2 cannot be both True, choose one of them")

if b_cublas:
    from datapack import cublas
    from basics import add_cudasync

    cw.start_if_def(f"_CUBLAS")
    cw.include("<cublas.h>")
    cw.end_if_def()

    cw.start_if_def(f"_CUBLAS")
    add_cudasync(cw)
    cw.end_if_def()

    for entry in cublas:
        add_this(cw, entry)

if b_cublas_v2:
    from datapack import cublas_v2
    from basics import ( add_cudasync,
                         add_cublas_handle_init,
                         add_cublas_handle_destroy,
                         add_cublas_types,
                       )

    cw.start_if_def(f"_CUBLAS")
    cw.include("<cublas_v2.h>")
    cw.end_if_def()

    cw.start_if_def(f"_CUBLAS")
    add_cublas_handle_init(cw)
    add_cublas_handle_destroy(cw)
    cw.end_if_def()

    cw.start_if_def(f"_CUBLAS")
    add_cudasync(cw)
    cw.end_if_def()

    cw.start_if_def(f"_CUBLAS")
    add_cublas_types(cw)
    cw.end_if_def()

    for entry in cublas_v2:
        add_this(cw, entry)

if b_cusolver:
    from datapack import cusolver
    from basics import ( add_cusolver_handle_init,
                         add_cusolver_handle_destroy,
                       )

    cw.start_if_def(f"_CUSOLVER")
    cw.include('"cusolverDn.h"')
    cw.end_if_def()

    cw.start_if_def(f"_CUSOLVER")
    add_cusolver_handle_init(cw)
    add_cusolver_handle_destroy(cw)
    cw.end_if_def()

    cw.start_if_def(f"_CUSOLVER")

    for entry in cusolver:
        add_this(cw, entry)

print(cw)
