import pyomo.environ as pyo
import math
import random

def solve_sphere_packing(N, r_s, R_c, H_c, timeout_seconds=120):
    """
    Attempts to find a feasible packing for N spheres using Pyomo and a non-linear solver.
    
    Args:
        N (int): The number of spheres to pack.
        r_s (float): Radius of the spheres.
        R_c (float): Radius of the cylinder.
        H_c (float): Height of the cylinder.
        timeout_seconds (int): Time limit for the solver.

    Returns:
        (bool, model): A tuple containing a boolean for feasibility and the solved model.
    """
    print(f"\n--- Checking for N = {N} spheres ---")
    model = pyo.ConcreteModel(f"SpherePacking_N{N}")

    # --- Sets ---
    model.I = pyo.RangeSet(1, N)
    model.J = pyo.RangeSet(1, N)

    # --- Parameters ---
    model.r_s = pyo.Param(initialize=r_s)
    model.R_c = pyo.Param(initialize=R_c)
    model.H_c = pyo.Param(initialize=H_c)
    
    # Radius of the circle containing the sphere centers
    center_R_max = model.R_c - model.r_s
    
    # --- Variables ---
    # Initialize variables with random values to help the solver find a starting point.
    model.x = pyo.Var(model.I, bounds=(-center_R_max, center_R_max), initialize=lambda m, i: random.uniform(-center_R_max, center_R_max))
    model.y = pyo.Var(model.I, bounds=(-center_R_max, center_R_max), initialize=lambda m, i: random.uniform(-center_R_max, center_R_max))
    model.z = pyo.Var(model.I, bounds=(model.r_s, model.H_c - model.r_s), initialize=lambda m, i: random.uniform(model.r_s, model.H_c - model.r_s))

    # --- Constraints ---
    # 1. Spheres must be inside the cylinder's radius
    def cylinder_radius_rule(m, i):
        return m.x[i]**2 + m.y[i]**2 <= center_R_max**2
    model.CylinderRadiusConstraint = pyo.Constraint(model.I, rule=cylinder_radius_rule)

    # 2. Spheres must not overlap
    def non_overlap_rule(m, i, j):
        if i < j:
            return (m.x[i] - m.x[j])**2 + (m.y[i] - m.y[j])**2 + (m.z[i] - m.z[j])**2 >= (2 * m.r_s)**2
        return pyo.Constraint.Skip
    model.NonOverlapConstraint = pyo.Constraint(model.I, model.J, rule=non_overlap_rule)

    # --- Objective ---
    # This is a feasibility problem, so the objective is trivial.
    model.obj = pyo.Objective(expr=1.0)

    # --- Solve ---
    # Use IPOPT solver. 'couenne' or 'scip' would be better for a guaranteed global solution.
    try:
        solver = pyo.SolverFactory('ipopt')
        solver.options['max_cpu_time'] = timeout_seconds
        results = solver.solve(model, tee=False) # Set tee=True to see detailed solver output

        if results.solver.termination_condition in (pyo.TerminationCondition.optimal, pyo.TerminationCondition.feasible):
            print(f"Result: Feasible solution found for N = {N}.")
            return True, model
        else:
            print(f"Result: Solver could not find a feasible solution for N = {N}. Status: {results.solver.termination_condition}")
            return False, None
    except Exception as e:
        print(f"An error occurred while solving for N = {N}: {e}")
        print("Please ensure you have Pyomo and a solver like IPOPT installed.")
        print("Try: pip install pyomo && conda install -c conda-forge ipopt")
        return False, None

if __name__ == '__main__':
    # Problem parameters
    sphere_radius = 2.0
    cylinder_radius = 10.0
    cylinder_height = 11.0

    # Start searching downwards from a known good estimate.
    # Literature suggests N=56 is the max, and N=57 is impossible.
    start_n = 56 
    
    found_solution = False
    for n_spheres in range(start_n, 0, -1):
        is_feasible, solution_model = solve_sphere_packing(
            N=n_spheres,
            r_s=sphere_radius,
            R_c=cylinder_radius,
            H_c=cylinder_height
        )
        if is_feasible:
            print("\n=======================================================")
            print(f"Maximum number of spheres found: {n_spheres}")
            print("=======================================================")
            print("\nFinal equation parameters and sphere coordinates:")
            print(f"Sphere radius (r_s): {pyo.value(solution_model.r_s)}")
            print(f"Cylinder radius (R_c): {pyo.value(solution_model.R_c)}")
            print(f"Cylinder height (H_c): {pyo.value(solution_model.H_c)}")
            print(f"Non-overlap distance squared >= (2*r_s)^2 = {(2*pyo.value(solution_model.r_s))**2}")
            print("-" * 55)
            for i in solution_model.I:
                x_val = pyo.value(solution_model.x[i])
                y_val = pyo.value(solution_model.y[i])
                z_val = pyo.value(solution_model.z[i])
                print(f"Sphere {i:>2}: x={x_val:8.4f}, y={y_val:8.4f}, z={z_val:8.4f}")
            print("=======================================================")
            found_solution = True
            break # Stop searching
    
    if not found_solution:
        print("\nCould not find a feasible solution within the tested range and time limits.")

<<<56>>>