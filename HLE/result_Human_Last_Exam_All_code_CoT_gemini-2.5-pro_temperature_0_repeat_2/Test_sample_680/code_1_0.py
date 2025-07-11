import math
from scipy.optimize import milp

def solve_tower_optimization():
    """
    Formulates and solves the tower placement optimization problem.
    """
    # --- Problem Definition ---
    # City dimensions
    city_x = 12  # km
    city_y = 11  # km
    total_area = city_x * city_y

    # Tower specifications
    # B1: {radius, cost}
    b1_spec = {'r': 1, 'cost': 1500}
    # B2: {radius, cost}
    b2_spec = {'r': 2, 'cost': 5000}

    # Budget
    budget = 45000

    # --- ILP Formulation ---
    # The decision variables are x = [b1, b2]
    # We want to maximize the coverage area: pi*r1^2*b1 + pi*r2^2*b2
    # Maximize: pi * (1*b1 + 4*b2)
    # Since pi is a constant, this is equivalent to maximizing: 1*b1 + 4*b2
    # scipy.optimize.milp performs minimization, so we minimize the negative of our objective.
    # Objective function coefficients (c):
    c = [-1, -4]

    # Constraints (A_ub @ x <= b_ub)
    # 1. Budget constraint: 1500*b1 + 5000*b2 <= 45000
    #    Simplified by dividing by 500: 3*b1 + 10*b2 <= 90
    # 2. Area constraint: pi*r1^2*b1 + pi*r2^2*b2 <= total_area
    #    pi*(1*b1 + 4*b2) <= 132
    #    1*b1 + 4*b2 <= 132 / pi  (approx 42.01)
    #    Since b1 and b2 are integers, 1*b1 + 4*b2 <= 42
    A_ub = [[3, 10], [1, 4]]
    b_ub = [90, 42]

    # Both variables must be integers
    integrality = [1, 1]

    # Bounds for variables (b1 >= 0, b2 >= 0)
    bounds = [(0, None), (0, None)]

    # --- Solve the problem ---
    res = milp(c=c, A_ub=A_ub, b_ub=b_ub, integrality=integrality, bounds=bounds)

    if not res.success:
        print("Optimization failed.")
        return

    # Extract results
    b1 = int(round(res.x[0]))
    b2 = int(round(res.x[1]))

    # --- Calculate Coverage Ratio ---
    covered_area = math.pi * (b1 * b1_spec['r']**2 + b2 * b2_spec['r']**2)
    coverage_ratio = (covered_area / total_area) * 100
    c = int(round(coverage_ratio))

    # --- Output Results ---
    # As requested, showing the numbers in the final equation for the coverage ratio
    print(f"Optimal solution: b1={b1}, b2={b2}")
    print(f"Coverage ratio equation: c = round((pi * ({b1} * {b1_spec['r']**2} + {b2} * {b2_spec['r']**2})) / {total_area} * 100)")
    print(f"Resulting coverage ratio (c): {c}%")
    
    # Final answer in the required format
    print("\nFormatted Answer:")
    print(f"{b1};{b2};{c}")


if __name__ == '__main__':
    solve_tower_optimization()
