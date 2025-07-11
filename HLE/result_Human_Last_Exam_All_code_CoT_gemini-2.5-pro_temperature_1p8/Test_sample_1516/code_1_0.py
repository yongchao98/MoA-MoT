import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design problem.
    """
    # Problem parameters
    total_members = 791
    sections = 61
    row_depth = 1.5  # meters
    initial_radius = 3.0 # meters

    # Step 1: Calculate the number of rows per section
    rows_per_section = total_members // sections

    print("--- Step 1: Seating Arrangement ---")
    print(f"The parliament has {total_members} members and {sections} sections.")
    print(f"For symmetry, each section will have {total_members} / {sections} = {rows_per_section} rows.")
    print("-" * 40)

    # Step 2: Formulate the visibility constraint inequality
    # The derivation shows that for visibility to be guaranteed for all rows, K must satisfy:
    # K < 4.5 * (j-1) * (j-2)
    # This must hold for j from 3 to 13 (the last row).
    print("--- Step 2: Visibility Constraint ---")
    print("The geometric analysis leads to the inequality for K:")
    print("K < 4.5 * (j-1) * (j-2)")
    print("This must hold for j=3 (the first possible obstruction) up to j=13.")
    print("-" * 40)

    # Step 3: Find the most restrictive case
    # The expression 4.5 * (j-1) * (j-2) is smallest when j is smallest.
    j_critical = 3
    print("--- Step 3: Finding the Tightest Constraint ---")
    print("The constraint is tightest for the smallest value of j, which is j = 3.")
    print(f"This represents the view from row 3 over the person in row 2.")
    print("-" * 40)
    
    # Step 4: Calculate the final value of K
    j_minus_1 = j_critical - 1
    j_minus_2 = j_critical - 2
    factor = 4.5
    
    limit_for_K = factor * j_minus_1 * j_minus_2

    print("--- Step 4: Final Calculation ---")
    print(f"We plug j = {j_critical} into the inequality:")
    print(f"K < {factor} * ({j_critical} - 1) * ({j_critical} - 2)")
    print(f"K < {factor} * {j_minus_1} * {j_minus_2}")
    print(f"K < {limit_for_K}")
    print("\nSince K must be an integer and strictly less than the calculated limit...")
    
    # Since K must be strictly less than 9, the largest integer it can be is 8.
    max_K = math.floor(limit_for_K - 1e-9)

    print(f"The maximum integer value for K is {max_K}.")
    print(f"\n<<<8>>>")

solve_parliament_design()