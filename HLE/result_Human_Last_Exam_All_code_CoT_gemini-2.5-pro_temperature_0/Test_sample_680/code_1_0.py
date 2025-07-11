import pulp
import math

def solve_wifi_optimization():
    """
    This function formulates and solves the WiFi tower placement optimization problem.
    """
    # 1. Define the problem: We want to maximize the coverage.
    prob = pulp.LpProblem("Maximize_Wifi_Coverage", pulp.LpMaximize)

    # 2. Define the decision variables: number of towers of each type.
    # b1: number of B1 towers (radius 1km, cost 1500)
    b1 = pulp.LpVariable("b1", lowBound=0, cat='Integer')
    # b2: number of B2 towers (radius 2km, cost 5000)
    b2 = pulp.LpVariable("b2", lowBound=0, cat='Integer')

    # 3. Define the objective function.
    # Maximizing total area pi*(1^2*b1 + 2^2*b2) is equivalent to maximizing (b1 + 4*b2).
    prob += b1 + 4 * b2, "Total Coverage Factor"

    # 4. Define the constraints.
    # Constraint 1: Budget constraint
    # 1500*b1 + 5000*b2 <= 45000
    prob += 1500 * b1 + 5000 * b2 <= 45000, "BudgetConstraint"

    # Constraint 2: Area constraint (non-overlapping towers)
    # The total area of circles must be less than or equal to the city's area.
    city_area = 12 * 11
    # pi*1^2*b1 + pi*2^2*b2 <= city_area
    prob += math.pi * (1**2 * b1 + 2**2 * b2) <= city_area, "AreaConstraint"

    # 5. Solve the problem.
    # The status of the solution is printed to the console (e.g., "Optimal").
    prob.solve()

    # 6. Extract the optimal values for b1 and b2.
    optimal_b1 = int(pulp.value(b1))
    optimal_b2 = int(pulp.value(b2))

    # 7. Calculate the coverage and ratio.
    # The final equation for coverage area is based on the optimal number of towers.
    coverage_area_b1 = math.pi * (1**2) * optimal_b1
    coverage_area_b2 = math.pi * (2**2) * optimal_b2
    total_coverage_area = coverage_area_b1 + coverage_area_b2
    
    # The final equation for the coverage ratio.
    coverage_ratio = (total_coverage_area / city_area) * 100
    
    # Round the coverage ratio to the nearest percentage.
    coverage_percentage = round(coverage_ratio)

    # 8. Print the final result in the specified format.
    # The final equation is represented by its results: b1, b2, and the final coverage percentage c.
    print(f"{optimal_b1};{optimal_b2};{coverage_percentage}")

if __name__ == "__main__":
    solve_wifi_optimization()