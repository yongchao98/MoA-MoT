import pyomo.environ as pyo
from pyomo.opt import SolverFactory, TerminationCondition
import random
import math
import sys

def solve_packing(N, R, H, r):
    """
    Attempts to solve the sphere packing feasibility problem for N spheres.
    Returns True if a feasible solution is found, False otherwise.
    """
    model = pyo.ConcreteModel(f"Sphere_Packing_{N}")

    # --- Parameters ---
    model.R = R
    model.H = H
    model.r = r
    model.N = N
    model.indices = pyo.RangeSet(1, model.N)
    
    # --- Variables ---
    # Initialize variables with random values within the allowed space
    # This helps the non-linear solver find a solution.
    def x_init(m, i):
        return random.uniform(-(m.R - m.r), m.R - m.r)
    def y_init(m, i):
        # Ensure x^2 + y^2 <= (R-r)^2
        max_y = math.sqrt((m.R - m.r)**2 - m.x[i]**2)
        return random.uniform(-max_y, max_y)
    def z_init(m, i):
        return random.uniform(m.r, m.H - m.r)

    model.x = pyo.Var(model.indices, initialize=x_init, bounds=(-(R-r), R-r))
    model.y = pyo.Var(model.indices, initialize=y_init, bounds=(-(R-r), R-r))
    model.z = pyo.Var(model.indices, initialize=z_init, bounds=(r, H-r))

    # --- Constraints ---
    # 1. Cylinder containment constraint (radial)
    def containment_rule(m, i):
        return m.x[i]**2 + m.y[i]**2 <= (m.R - m.r)**2
    model.containment = pyo.Constraint(model.indices, rule=containment_rule)

    # 2. Non-overlapping constraint for each pair of spheres
    # To avoid duplicate constraints (i,j) vs (j,i), we define the rule over a set of pairs.
    model.pairs = pyo.Set(initialize=[(i,j) for i in model.indices for j in model.indices if i < j])
    def non_overlap_rule(m, i, j):
        return (m.x[i] - m.x[j])**2 + (m.y[i] - m.y[j])**2 + (m.z[i] - m.z[j])**2 >= (2 * m.r)**2
    model.non_overlap = pyo.Constraint(model.pairs, rule=non_overlap_rule)

    # --- Objective ---
    # This is a feasibility problem, so we can set a constant objective.
    model.objective = pyo.Objective(expr=0)

    # --- Solve ---
    try:
        # We need a non-linear solver like 'ipopt'
        solver = SolverFactory('ipopt')
        results = solver.solve(model, tee=False) # tee=True to see solver output
    except ApplicationError:
        print("IPOPT solver not found. Please make sure it is installed and in your system's PATH.")
        sys.exit(1)


    # Check solver termination condition
    if results.solver.termination_condition in [TerminationCondition.optimal, TerminationCondition.locallyOptimal]:
        # A feasible solution was found
        # We store the solution coordinates in the model object itself for later printing
        model.solution_found = True
        model.solution_coords = {i: (pyo.value(model.x[i]), pyo.value(model.y[i]), pyo.value(model.z[i])) for i in model.indices}
        return True, model
    else:
        # No feasible solution found by the solver
        model.solution_found = False
        return False, None


if __name__ == '__main__':
    # Problem parameters
    R_cyl = 10  # Cylinder radius
    H_cyl = 11  # Cylinder height
    r_sph = 2   # Sphere radius

    # Use binary search to find the maximum N
    low = 1
    high = 100 # A reasonable upper bound based on volume
    best_n = 0
    best_model = None

    print(f"Searching for the maximum number of spheres (radius={r_sph}) that can be packed in a cylinder (radius={R_cyl}, height={H_cyl}).")
    print("This may take a few minutes...")

    while low <= high:
        mid = (low + high) // 2
        if mid == 0: break
        print(f"Testing N = {mid}... ", end="")
        # Run solver multiple times with random restarts to improve chance of finding a solution
        is_feasible = False
        solved_model = None
        for attempt in range(3): # 3 attempts for robustness
            feasible, model_instance = solve_packing(N=mid, R=R_cyl, H=H_cyl, r=r_sph)
            if feasible:
                is_feasible = True
                solved_model = model_instance
                break
        
        if is_feasible:
            print("Feasible solution found.")
            best_n = mid
            best_model = solved_model
            low = mid + 1  # Try for more spheres
        else:
            print("No feasible solution found.")
            high = mid - 1 # Try for fewer spheres

    print("\n" + "="*40)
    print("Search complete.")
    print(f"Cylinder Radius (R): {R_cyl}")
    print(f"Cylinder Height (H): {H_cyl}")
    print(f"Sphere Radius (r): {r_sph}")
    print(f"\nThe maximum number of spheres that can be packed is: {best_n}")
    print("="*40)

    if best_model and best_model.solution_found:
        print("\nCoordinates of the sphere centers for the optimal packing found:")
        for i, coords in sorted(best_model.solution_coords.items()):
            print(f"Sphere {i:>2}: (x={coords[0]:>8.4f}, y={coords[1]:>8.4f}, z={coords[2]:>8.4f})")
    
    print(f"\n<<<{best_n}>>>")
