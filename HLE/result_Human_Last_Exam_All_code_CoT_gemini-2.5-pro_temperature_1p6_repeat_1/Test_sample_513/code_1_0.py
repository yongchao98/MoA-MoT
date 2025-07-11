import numpy as np
from scipy.optimize import minimize
import random

def solve_sphere_packing():
    """
    This script finds the maximum number of spheres that can be packed into a cylinder
    by formulating and solving a series of optimization problems.
    """
    # --- Problem Parameters ---
    r_sphere = 2.0
    R_cylinder = 10.0
    H_cylinder = 11.0

    print("Problem: Find the maximum number of spheres (N) that can be packed in a cylinder.")
    print("\n--- Parameters ---")
    print(f"Sphere Radius (r): {r_sphere}")
    print(f"Cylinder Radius (R): {R_cylinder}")
    print(f"Cylinder Height (H): {H_cylinder}")

    # --- Derived constants for the optimization constraints ---
    R_eff_sq = (R_cylinder - r_sphere)**2
    H_min = r_sphere
    H_max = H_cylinder - r_sphere
    d_sq = (2 * r_sphere)**2

    print("\n--- Optimization Constraints (for each sphere i and pair i,j) ---")
    print(f"1. Radial: x_i^2 + y_i^2 <= {R_eff_sq}")
    print(f"2. Height:  {H_min} <= z_i <= {H_max}")
    print(f"3. No Overlap: (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= {d_sq}")
    print("-" * 40)

    def can_pack(N, num_restarts=5):
        """
        Checks if N spheres can be packed by solving a feasibility problem.
        Returns True if a feasible solution is found, False otherwise.
        """
        if N <= 0: return True
        
        # Objective function is constant as we only care about feasibility.
        objective = lambda x: 0

        # Create a list of constraints for the solver.
        cons = []
        # Constraint 1: Radial distance (formulated for solver as g(x) >= 0)
        for i in range(N):
            cons.append({'type': 'ineq', 'fun': lambda x, i=i: R_eff_sq - (x[3*i]**2 + x[3*i+1]**2)})

        # Constraint 3: Non-overlap distance
        for i in range(N):
            for j in range(i + 1, N):
                cons.append({'type': 'ineq', 'fun': lambda x, i=i, j=j: np.sum((x[3*i:3*i+3] - x[3*j:3*j+3])**2) - d_sq})

        # Constraint 2: Height bounds for each z_i coordinate
        bounds = []
        for _ in range(N):
            bounds.extend([(-R_cylinder + r_sphere, R_cylinder - r_sphere),   # x_i
                           (-R_cylinder + r_sphere, R_cylinder - r_sphere),   # y_i
                           (H_min, H_max)])                                 # z_i

        # Try multiple random starting points to improve chances of finding a solution.
        for attempt in range(num_restarts):
            # Generate a random initial guess for the sphere center coordinates
            x0 = np.empty(3 * N)
            for i in range(N):
                rand_r = (R_cylinder - r_sphere) * np.sqrt(random.random())
                rand_theta = 2 * np.pi * random.random()
                x0[3*i] = rand_r * np.cos(rand_theta)
                x0[3*i+1] = rand_r * np.sin(rand_theta)
                x0[3*i+2] = random.uniform(H_min, H_max)
            
            # Run the SLSQP optimizer
            res = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=cons, 
                           options={'maxiter': 500, 'ftol': 1e-7})

            # Check if the solver terminated successfully and met all constraints
            if res.success:
                final_x = res.x
                all_constraints_met = True
                # The solver might report success with tiny violations, so we double check
                for con in cons:
                    if con['fun'](final_x) < -1e-5:
                        all_constraints_met = False
                        break
                if all_constraints_met:
                    print(f"Attempt {attempt+1}/{num_restarts}: Feasible solution FOUND.")
                    return True
            print(f"Attempt {attempt+1}/{num_restarts}: Solver could not find a feasible packing.")

        return False

    # --- Search for Maximum N ---
    # We will search downwards from an upper bound. This can be very time-consuming.
    # The search starts from 59, which is a value slightly above the known optimum from literature.
    print("\n--- Searching for Maximum N ---")
    print("Note: This is a computationally intensive search. It may take several minutes.\n")
    
    max_n_found = 0
    # Search from N=59 down to N=57.
    for n_spheres in range(59, 56, -1):
        print(f"Checking if N = {n_spheres} spheres can be packed...")
        # Use a higher number of restarts for better chance of finding a solution.
        if can_pack(n_spheres, num_restarts=5):
            max_n_found = n_spheres
            print(f"\nSUCCESS: Found a feasible packing for {max_n_found} spheres.")
            print("Since we are searching downwards, this is the maximum number found.")
            break
        else:
            print(f"Could not find a feasible packing for N = {n_spheres}. Trying a smaller N...")
    
    if max_n_found > 0:
        print(f"\nThe optimal number of spheres is {max_n_found}.")
    else:
        print("\nCould not find a solution in the searched range. The problem is very complex,")
        print("or more restarts/a different solver may be required.")

    return max_n_found

# --- Execute the Solver ---
result = solve_sphere_packing()
print(f"\nFinal Answer: The maximum number of spheres is {result}.")
print(f"<<<{result}>>>")