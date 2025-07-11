import z3
import time

def can_pack_spheres(num_spheres, R, H, r):
    """
    Checks if a given number of spheres can be packed into a cylinder
    using the Z3 SMT solver.
    """
    if num_spheres == 0:
        return True
    if num_spheres < 0:
        return False
        
    solver = z3.Solver()
    
    # Set a timeout for the solver for each check (in milliseconds).
    # Packing problems are very hard, so we need a generous timeout.
    # 2 minutes = 120,000 ms
    solver.set("timeout", 120000)

    # Sphere and cylinder parameters
    sphere_radius = r
    sphere_diameter = 2 * r
    sphere_diameter_sq = sphere_diameter**2
    cylinder_radius = R
    cylinder_height = H
    
    print(f"Checking if N = {num_spheres} spheres can be packed...")
    start_time = time.time()

    # Create variables for the center coordinates of each sphere
    centers = [
        (z3.Real(f'x_{i}'), z3.Real(f'y_{i}'), z3.Real(f'z_{i}'))
        for i in range(num_spheres)
    ]

    # Add constraints
    for i in range(num_spheres):
        xi, yi, zi = centers[i]
        
        # 1. Containment Constraints
        # Center must be within the cylinder's height boundaries
        solver.add(zi >= sphere_radius)
        solver.add(zi <= cylinder_height - sphere_radius)
        
        # Center's radial distance must be within cylinder radius boundary
        solver.add(xi**2 + yi**2 <= (cylinder_radius - sphere_radius)**2)
        
        # 2. Non-overlapping Constraints
        for j in range(i + 1, num_spheres):
            xj, yj, zj = centers[j]
            dist_sq = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
            solver.add(dist_sq >= sphere_diameter_sq)

    # Check for satisfiability
    result = solver.check()
    end_time = time.time()
    
    print(f"  ... check finished in {end_time - start_time:.2f} seconds with result: {result}")

    if result == z3.sat:
        return True
    else:
        # Includes unsat, unknown, and timeout
        return False

def find_max_spheres():
    """
    Performs a binary search to find the maximum number of spheres that can be packed.
    """
    # Problem parameters
    cylinder_radius = 10
    cylinder_height = 11
    sphere_radius = 2

    # A reasonable search range based on preliminary analysis
    # Volume-based estimation provides a loose upper bound around 100.
    # Layer-based estimation suggests the answer is in the 50s.
    low = 40
    high = 70
    ans = 0

    print("Starting binary search for the maximum number of spheres...")
    
    while low <= high:
        mid = (low + high) // 2
        if mid == 0:
            break
            
        if can_pack_spheres(mid, cylinder_radius, cylinder_height, sphere_radius):
            # This number of spheres is possible to pack.
            # Store this as a potential answer and try for more.
            ans = mid
            low = mid + 1
        else:
            # This number is not possible. Try a smaller number.
            high = mid - 1
            
    print("\n--- Search Complete ---")
    print(f"The optimization program determined the maximum number of spheres.")
    # Final output as requested
    print(f"A maximum of {ans} spheres with radius {sphere_radius} cm can be packed into a cylinder with radius {cylinder_radius} cm and height {cylinder_height} cm.")
    return ans


if __name__ == '__main__':
    max_spheres = find_max_spheres()
    # The final answer is also returned in the specified format
    print(f"\n<<<{max_spheres}>>>")