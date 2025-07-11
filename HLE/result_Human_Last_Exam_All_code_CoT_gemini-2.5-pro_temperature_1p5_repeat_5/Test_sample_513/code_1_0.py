import math
import random
from pyomo.environ import ConcreteModel, RangeSet, Var, Constraint, Objective, SolverFactory, SolverStatus, TerminationCondition, NonNegativeReals

def solve_sphere_packing():
    """
    Solves the sphere packing problem for the given dimensions
    using a binary search over the number of spheres.
    """
    R_c = 10  # Cylinder radius
    H_c = 11  # Cylinder height
    r_s = 2   # Sphere radius

    # Store the model from the last successful packing to retrieve coordinates
    last_successful_model = None

    def can_pack(N):
        """
        Checks if N spheres can be packed by solving a feasibility problem.
        Returns True if feasible, False otherwise.
        """
        nonlocal last_successful_model
        model = ConcreteModel(f"Sphere_Packing_N_{N}")

        # --- Sets ---
        model.I = RangeSet(1, N)
        model.J = RangeSet(1, N)

        # --- Parameters ---
        model.R_c = R_c
        model.H_c = H_c
        model.r_s = r_s

        # --- Variables ---
        # Initialize with random values to help the solver find a starting point
        def get_rand_x(model): return random.uniform(-(R_c - r_s), R_c - r_s)
        def get_rand_y(model): return random.uniform(-(R_c - r_s), R_c - r_s)
        def get_rand_z(model): return random.uniform(r_s, H_c - r_s)

        model.x = Var(model.I, initialize=get_rand_x, bounds=(-(R_c - r_s), R_c - r_s))
        model.y = Var(model.I, initialize=get_rand_y, bounds=(-(R_c - r_s), R_c - r_s))
        model.z = Var(model.I, initialize=get_rand_z, bounds=(r_s, H_c - r_s))

        # --- Constraints ---
        # Constraint 1: All spheres must be within the cylinder's radial boundary
        def cylinder_rule(model, i):
            return model.x[i]**2 + model.y[i]**2 <= (model.R_c - model.r_s)**2
        model.CylinderConstraint = Constraint(model.I, rule=cylinder_rule)

        # Constraint 2: Non-overlapping spheres
        def non_overlap_rule(model, i, j):
            if i >= j:
                return Constraint.Skip
            return (model.x[i] - model.x[j])**2 + (model.y[i] - model.y[j])**2 + (model.z[i] - model.z[j])**2 >= (2 * model.r_s)**2
        model.NonOverlapConstraint = Constraint(model.I, model.J, rule=non_overlap_rule)

        # --- Objective ---
        # We only care about feasibility, so the objective is a dummy one.
        model.obj = Objective(expr=1.0)

        # --- Solve ---
        # Using 'ipopt' which is suitable for non-linear problems.
        solver = SolverFactory('ipopt')
        solver.options['tol'] = 1e-6 # Set solver tolerance
        # The 'tee=True' option shows solver output, which can be useful for debugging.
        # Set to False for a cleaner output.
        results = solver.solve(model, tee=False)

        # --- Check Result ---
        if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition in [TerminationCondition.optimal, TerminationCondition.feasible]):
            # If a feasible solution is found, store the model
            last_successful_model = model.clone()
            return True
        return False

    # --- Binary Search for the maximum number of spheres ---
    # The absolute upper bound is determined by volume ratio, which is around 103.
    low = 1
    high = 104
    max_n = 0

    print("Starting binary search to find the maximum number of spheres...")
    print("This may take several minutes depending on your system's performance.")
    while low <= high:
        mid = (low + high) // 2
        if mid == 0:
            break
        print(f"Testing feasibility for N = {mid} spheres... ", end="")
        if can_pack(mid):
            print(f"SUCCESSFUL.")
            max_n = mid  # Found a feasible solution, try for more
            low = mid + 1
        else:
            print(f"FAILED.")
            high = mid - 1 # Could not find a feasible solution, try for less

    print("\n" + "="*50)
    print("                OPTIMIZATION COMPLETE")
    print("="*50)

    # --- Print Final Results ---
    if max_n > 0:
        print(f"The maximum number of spheres that can be packed is {max_n}.")
        vol_spheres = max_n * (4/3) * math.pi * r_s**3
        vol_cylinder = math.pi * R_c**2 * H_c
        packing_density = vol_spheres / vol_cylinder
        print(f"This corresponds to a packing density of {packing_density:.4f}.")
        print("-" * 50)
        print("Final Equation Solution:")
        print(f"The problem is to Maximize N, with the final result N = {max_n}.")
        print("The solution provides the coordinates for each sphere center such that")
        print("all constraints (non-overlap and containment) are satisfied.")
        print("\nSolved sphere center coordinates:")
        
        if last_successful_model is not None:
            for i in sorted(last_successful_model.I):
                print(f"  Sphere {i:2d}: (x={last_successful_model.x[i].value:8.4f}, y={last_successful_model.y[i].value:8.4f}, z={last_successful_model.z[i].value:8.4f})")
    else:
        print("Could not find a feasible packing for any number of spheres.")

    return max_n

# --- Run the solver ---
if __name__ == "__main__":
    final_answer = solve_sphere_packing()
    print(f"\n<<< {final_answer} >>>")
