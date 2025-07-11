import pyomo.environ as pyo
import random
import sys
import math

def solve_sphere_packing():
    """
    This script attempts to solve the sphere packing problem for a given cylinder and sphere size.
    It formulates the problem as a non-linear feasibility program and uses the Pyomo modeling
    language with the IPOPT solver to find a solution.

    The script iterates downwards from a high estimate for the number of spheres (N) and stops
    at the first N for which it can find a valid packing arrangement. Due to the non-convex
    nature of the problem, the solver is not guaranteed to find a solution even if one exists.
    To mitigate this, each N is attempted multiple times with different random initial positions.
    """

    # --- Problem Parameters ---
    R_cylinder = 10.0  # Cylinder radius
    H_cylinder = 11.0  # Cylinder height
    r_sphere = 2.0     # Sphere radius

    # Effective dimensions for sphere centers
    R_eff = R_cylinder - r_sphere
    Z_min = r_sphere
    Z_max = H_cylinder - r_sphere

    print("--- Problem Formulation ---")
    print(f"Objective: Find the maximum number of spheres (N) that can be packed into a cylinder.")
    print("\nParameters:")
    print(f"  - Sphere Radius (r): {r_sphere}")
    print(f"  - Cylinder Radius (R): {R_cylinder}")
    print(f"  - Cylinder Height (H): {H_cylinder}")
    print("\nVariables:")
    print(f"  - N: Number of spheres (integer, to be maximized)")
    print(f"  - (x_i, y_i, z_i): Center coordinates of sphere i, for i = 1 to N")
    print("\nConstraints:")
    print(f"  1. Cylinder Containment: x_i^2 + y_i^2 <= (R - r)^2")
    print(f"                           r <= z_i <= H - r")
    print(f"  2. Non-Overlap: (x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2 >= (2*r)^2 for i != j")
    print("-------------------------\n")
    
    # --- Solver Configuration ---
    # Search for N by starting high and going down. Based on heuristics, the answer is likely in the 55-60 range.
    start_n = 60
    end_n = 50
    num_trials_per_n = 5 # Number of random initializations to try for each N

    # Check if a solver is available
    try:
        solver = pyo.SolverFactory('ipopt')
        if not solver.available():
             raise OSError("IPOPT solver not found.")
    except (OSError, pyo.common.errors.ApplicationError):
        print("ERROR: The 'ipopt' solver is required but not found in the system's PATH.")
        print("Please install IPOPT. For example, on Ubuntu/Debian: sudo apt-get install coinor-ipopt")
        print("Or via conda: conda install -c conda-forge ipopt")
        return
        
    best_n_found = 0
    solution_model = None

    for n in range(start_n, end_n - 1, -1):
        print(f"Trying to pack N = {n} spheres...")
        is_packable = False
        for trial in range(num_trials_per_n):
            model = pyo.ConcreteModel(f"SpherePacking_N{n}")
            
            # --- Model Sets ---
            model.I = pyo.RangeSet(1, n)
            # Create a set of unique pairs (i, j) where i < j
            model.P = pyo.Set(initialize=[(i, j) for i in model.I for j in model.I if i < j])

            # --- Model Variables ---
            model.x = pyo.Var(model.I, bounds=(-R_eff, R_eff))
            model.y = pyo.Var(model.I, bounds=(-R_eff, R_eff))
            model.z = pyo.Var(model.I, bounds=(Z_min, Z_max))
            
            # --- Constraints ---
            # 1. Cylinder containment
            @model.Constraint(model.I)
            def cylinder_constraint(m, i):
                return m.x[i]**2 + m.y[i]**2 <= R_eff**2

            # 2. Non-overlapping spheres
            @model.Constraint(model.P)
            def non_overlap_constraint(m, i, j):
                return (m.x[i] - m.x[j])**2 + (m.y[i] - m.y[j])**2 + (m.z[i] - m.z[j])**2 >= (2 * r_sphere)**2

            # --- Objective Function ---
            # This is a feasibility problem, so we provide a dummy objective.
            model.obj = pyo.Objective(expr=1.0)
            
            # --- Initialization ---
            # Provide a random starting point for the solver.
            for i in model.I:
                angle = random.uniform(0, 2 * math.pi)
                radius = random.uniform(0, R_eff)
                model.x[i] = radius * math.cos(angle)
                model.y[i] = radius * math.sin(angle)
                model.z[i] = random.uniform(Z_min, Z_max)

            # --- Solve ---
            try:
                results = solver.solve(model, tee=False)
                # Check solver status
                if results.solver.termination_condition in [pyo.TerminationCondition.optimal, pyo.TerminationCondition.feasible]:
                    print(f"  Trial {trial + 1}/{num_trials_per_n}: Success! Found a valid packing for N = {n}.")
                    best_n_found = n
                    solution_model = model.clone()
                    is_packable = True
                    break # Stop trials for this N, as we have found a success.
                else:
                    sys.stdout.write(f"\r  Trial {trial + 1}/{num_trials_per_n}: Solver could not find a solution.")
                    sys.stdout.flush()

            except Exception as e:
                print(f"  Trial {trial + 1}/{num_trials_per_n}: An error occurred during solving: {e}")
        
        print() # Newline after trials for this N
        if is_packable:
            # We found the highest N in our search range that is packable, so we can stop.
            break

    # --- Print Final Answer ---
    if best_n_found > 0:
        print("\n--- SOLUTION ---")
        print(f"The maximum number of spheres found is: {best_n_found}")
        print("\nFinal Equation Parameters:")
        print(f"Number of Spheres (N): {best_n_found}")
        print(f"Sphere Radius (r): {r_sphere}")
        print(f"Cylinder Radius (R): {R_cylinder}")
        print(f"Cylinder Height (H): {H_cylinder}")
        print("\nFound sphere center coordinates (x, y, z):")
        for i in solution_model.I:
            print(f"  Sphere {i:2d}: ({solution_model.x[i].value:8.4f}, {solution_model.y[i].value:8.4f}, {solution_model.z[i].value:8.4f})")
        
        print(f"\n<<<{best_n_found}>>>")
    else:
        print(f"\n--- NO SOLUTION FOUND ---")
        print(f"Could not find a valid packing for N in the range [{start_n}, {end_n}].")
        print("You could try expanding the search range or increasing the number of trials.")
        print(f"\n<<<'Not found'>>>")

if __name__ == '__main__':
    solve_sphere_packing()