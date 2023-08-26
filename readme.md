# Hypersphere
This crate implements some simple rotation and projection primitives for 4D geometry on top of `glam`.

## Projection
Implements 4D to 3D projection on the surface of a sphere by way of stereographic projection.

It allows you to project _both_ points and vectors tangent to the sphere at a point through a stereographic projection.
The latter is useful when embedding 3D geometry on the surface of a hypersphere and ensuring that normal vectors remain
normal vectors under projection (recall that stereographic projection is angle-preserving).

## 4D Rotations
Implements a double-quaternion representation of 4D rotations. 

Includes:
- Rotation through basis planes (`XY`, `XZ`, etc.).
- Rotation through arbitrary pairs of planes specified by orthonormal vectors.
- Minimal rotations from one point to another.
- Cayley's decomposition of arbitrary 4D rotation matrices into this crate's representation.
- Slerp, inherited from quaternions.

## Basis Utilities
Includes functions to:
- Construct an arbitrary orthogonal vector to another vector.
- Construct an arbitrary orthogonal vector to two vectors.
- Construct a scaled version of the orthogonal vector to three vectors.
- Construct an orthonormal basis given two vectors that span a plane. 

## Sample Data
Implements a simple algorithm to generate the 600-cell's vertices (not indices), as useful sample data.
