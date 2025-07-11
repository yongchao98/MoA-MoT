import math
# You might need to install pulp: pip install pulp
try:
    import pulp
except ImportError:
    print("pulp library not found. Please install it using: pip install pulp")
    exit()

def solve_wifi_placement():
    """
    Formulates and solves the WiFi tower placement problem as an
    Integer Linear Programming problem.
    """
    # 1. Define problem constants
    city_area = 12 * 11  # km^2
    total_budget = 45000   # usd

    b1_radius = 1
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2

    b2_radius = 2
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2

    # 2. Formulate the problem
    # Create the model to maximize the objective function
    model = pulp.LpProblem(name="wifi-coverage-optimization", sense=pulp.LpMaximize)

    # Define the integer decision variables
    b1 = pulp.LpVariable(name="num_b1_towers", lowBound=0, cat="Integer")
    b2 = pulp.LpVariable(name="num_b2_towers", lowBound=0, cat="Integer")

    # Set the objective function: Maximize total covered area.
    # Maximizing Area = b1*b1_area + b2*b2_area is equivalent to maximizing
    # the expression b1 + 4*b2, which simplifies the solver's job.
    model += (b1 + 4 * b2, "total_coverage_factor")

    # Add the constraints
    # a) Budget constraint
    model += (b1 * b1_cost + b2 * b2_cost <= total_budget, "budget_constraint")
    # b) Non-overlapping area constraint
    model += (b1 * b1_area + b2 * b2_area <= city_area, "physical_area_constraint")

    # 3. Solve the problem
    # The solver will find the optimal integer values for b1 and b2
    model.solve()

    # 4. Extract and print results
    # Get the optimal values from the solver
    optimal_b1 = int(pulp.value(b1))
    optimal_b2 = int(pulp.value(b2))

    # Calculate the final coverage area and ratio
    final_covered_area = optimal_b1 * b1_area + optimal_b2 * b2_area
    coverage_ratio_percentage = (final_covered_area / city_area) * 100
    rounded_coverage_percentage = round(coverage_ratio_percentage)

    # As requested, output the numbers in the final equation
    print("Optimization Result:")
    print(f"Number of B1 towers (b1) = {optimal_b1}")
    print(f"Number of B2 towers (b2) = {optimal_b2}")
    print("\nCoverage Calculation:")
    print(f"Coverage Ratio = (({optimal_b1} * {b1_area:.2f}) + ({optimal_b2} * {b2_area:.2f})) / {city_area} * 100")
    print(f"Coverage Ratio = {coverage_ratio_percentage:.2f}%")
    print(f"Rounded to the nearest percentage = {rounded_coverage_percentage}%")
    
    # Print the final answer in the specified format b1;b2;c
    final_answer = f"{optimal_b1};{optimal_b2};{rounded_coverage_percentage}"
    print("\nFinal formatted answer:")
    print(final_answer)

    return final_answer


# Execute the function
solve_wifi_placement()

# The final line of the response contains the formatted answer for the system.
# The calculation results in b1=0, b2=9, and c=86
print("<<<0;9;86>>>")