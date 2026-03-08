from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from skimage.measure import marching_cubes as sk_marching_cubes


@dataclass
class IsosurfaceMesh:
    """Triangle mesh from isosurface extraction."""
    vertices: np.ndarray
    normals: np.ndarray
    faces: np.ndarray
    signs: np.ndarray


def marching_cubes(volume: np.ndarray, isovalue: float,
                   grid_origin: tuple[float, float, float],
                   grid_spacing: tuple[float, float, float]) -> IsosurfaceMesh:
    """Extract an isosurface from a 3D scalar field using marching cubes.

    Algorithm overview:
      1. For each voxel (2x2x2 cube of grid points), classify each corner
         as inside (value > isovalue) or outside.
      2. The 8-bit classification index selects a triangulation from a
         precomputed lookup table of 256 entries (15 unique cases).
      3. Interpolate vertex positions along edges where the isosurface
         crosses (linear interpolation between corner values).
      4. Compute vertex normals from the gradient of the scalar field.

    Uses sk_marching_cubes(volume, level=isovalue, spacing=grid_spacing)
    which returns (verts, faces, normals, values).

    After extraction, translate vertices to world coordinates by adding
    np.array(grid_origin). Set all signs to +1.0 (float32).

    Cast vertices and normals to float32, faces to uint32.

    Gotcha: Fails if isovalue is outside the volume's value range.
    Verify: volume.min() < isovalue < volume.max().

    Args:
        volume: 3D scalar field, shape (nx, ny, nz).
        isovalue: Scalar threshold for the isosurface.
        grid_origin: (x0, y0, z0) world coordinates of the grid corner.
        grid_spacing: (dx, dy, dz) spacing between grid points.

    Returns:
        IsosurfaceMesh with float32 vertices/normals, uint32 faces,
        and float32 signs (all +1.0).
    """
    ...


def extract_dual_isosurface(volume: np.ndarray, isovalue: float,
                            grid_origin: tuple[float, float, float],
                            grid_spacing: tuple[float, float, float]
                            ) -> IsosurfaceMesh:
    """Extract both +isovalue and -isovalue isosurfaces (dual lobes).

    Calls marching_cubes twice:
      1. Positive lobe: marching_cubes(volume, +isovalue, ...) -> signs = +1.
      2. Negative lobe: marching_cubes(volume, -isovalue, ...) -> set signs = -1.

    Combines meshes:
      - Concatenate vertices, normals, signs from both lobes.
      - Offset negative mesh's face indices by n_verts_positive.
      - Concatenate face arrays.

    Gotcha: For negative isosurface, verify volume.min() < -isovalue.
    At saddle points or very small isovalues near zero, marching cubes
    can produce degenerate triangles.

    Args:
        volume: 3D scalar field, shape (nx, ny, nz).
        isovalue: Positive scalar threshold; surfaces at +isovalue and -isovalue.
        grid_origin: (x0, y0, z0) world coordinates of the grid corner.
        grid_spacing: (dx, dy, dz) spacing between grid points.

    Returns:
        Combined IsosurfaceMesh with both lobes, signs indicating +1/-1.
    """
    ...
