import math

def solve_horse_water_problem(n, m):
    """
    Calculates the maximum water left and prints the equation.

    Args:
        n: The initial number of 100-liter water loads.
        m: The total distance to the destination in kilometers.
    """
    if n * 100 < m:
        print("Error: Not enough water to travel the distance directly even with one trip.")
        return

    # Find p, the number of trips required for the final leg of the journey.
    # It's the stage k where the destination m is reached.
    p = 1  # Default value if m is very large
    cummulative_distance = 0.0
    for k in range(n, 1, -1):
        segment_distance = 100.0 / (2 * k - 1)
        if m <= cummulative_distance + segment_distance:
            p = k
            break
        cummulative_distance += segment_distance

    # Build the string for the summation part of the equation
    sum_terms_list = []
    for k_val in range(p + 1, n + 1):
        sum_terms_list.append(f"100/(2*{k_val}-1)")
    
    sum_expr = " + ".join(sum_terms_list)

    if not sum_expr:
        sum_expr = "0"
    # Add parentheses if there are multiple terms in the sum
    elif len(sum_terms_list) > 1:
        sum_expr = f"({sum_expr})"

    # Construct the final equation string with all numbers and the summation expanded
    final_equation = f"Maximum Water Left = {p}*100 - ({m} - {sum_expr}) * (2*{p}-1)"

    # Print the final equation as requested
    print(final_equation)
    
    # Calculate the final numeric answer
    sum_val = 0.0
    for k_val in range(p + 1, n + 1):
      sum_val += 100.0 / (2 * k_val - 1)

    water_left = p * 100 - (m - sum_val) * (2 * p - 1)

    # Output the numeric answer in the specified format
    print(f"\n<<<{water_left:.2f}>>>")

# --- User-defined variables ---
# n: represents n*100 liters of initial water
n = 3
# m: represents the distance in kilometers
m = 50
# -----------------------------

# Assume n and m are provided. For this example, let's use n=3 and m=50.
# You can change these values to solve for different scenarios.
print(f"Solving for n={n} (initial water = {n*100}L) and m={m}km distance.\n")
solve_horse_water_problem(n, m)