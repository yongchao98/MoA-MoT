import pulp
import math

def solve_wifi_optimization():
    """
    Formulates and solves the WiFi tower placement problem to maximize coverage.
    """
    # 1. Define problem constants
    city_area = 12 * 11  # Total area in km^2
    total_budget = 45000  # Total budget in USD

    # B1 tower specifications
    b1_radius = 1
    b1_area = math.pi * (b1_radius**2)
    b1_cost = 1500

    # B2 tower specifications
    b2_radius = 2
    b2_area = math.pi * (b2_radius**2)
    b2_cost = 5000

    # 2. Formulate the Integer Linear Programming problem
    # Create a maximization problem
    prob = pulp.LpProblem("Maximize_Wifi_Coverage", pulp.LpMaximize)

    # Define decision variables (must be non-negative integers)
    b1 = pulp.LpVariable("num_b1_towers", lowBound=0, cat='Integer')
    b2 = pulp.LpVariable("num_b2_towers", lowBound=0, cat='Integer')

    # Define the objective function: Maximize the total coverage area
    prob += b1_area * b1 + b2_area * b2, "Total_Coverage_Area"

    # Define the constraints
    # Constraint 1: The total cost must be within the budget
    prob += b1_cost * b1 + b2_cost * b2 <= total_budget, "Budget_Constraint"
    
    # Constraint 2: The sum of tower areas cannot exceed the city's total area
    # This is a necessary condition for non-overlapping placement.
    prob += b1_area * b1 + b2_area * b2 <= city_area, "Area_Constraint"

    # 3. Solve the problem
    # The solver will find the optimal integer values for b1 and b2
    prob.solve()

    # 4. Extract and present the results
    # Get the optimal number of each tower type
    optimal_b1 = int(b1.varValue)
    optimal_b2 = int(b2.varValue)

    # Calculate the resulting maximum coverage area
    max_coverage_area = (b1_area * optimal_b1) + (b2_area * optimal_b2)

    # Calculate the coverage ratio and round it to the nearest percentage
    coverage_ratio_percentage = round((max_coverage_area / city_area) * 100)

    # Print the final result in the format b1;b2;c
    print(f"{optimal_b1};{optimal_b2};{coverage_ratio_percentage}")

if __name__ == "__main__":
    # To run this code, you might need to install the pulp library:
    # pip install pulp
    solve_wifi_optimization()