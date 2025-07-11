import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left and expresses the answer as an equation.

    Args:
        n (int): The number of 100-liter barrels of water.
        m (int): The total distance to travel in kilometers.
    """

    if n * 100 < m:
        print("Error: Not enough water to reach the destination.")
        return

    # Step 1: Find k, the number of full 100-liter barrels consumed to create depots.
    k = 0
    dist_covered = 0.0
    for i in range(1, n):
        # Number of barrels to move in this stage
        barrels_to_move = n - i + 1
        # Number of trips (forward and back) required
        num_trips = 2 * barrels_to_move - 1
        # Distance covered by consuming one full barrel (100 liters)
        dist_for_one_barrel = 100.0 / num_trips

        if dist_covered + dist_for_one_barrel <= m:
            dist_covered += dist_for_one_barrel
            k = i
        else:
            # The remaining distance is less than what would be covered by a full barrel,
            # so this is the final leg of the journey.
            break

    # Step 2: Construct the string for the summation part of the equation.
    sum_terms_list = []
    if k > 0:
        for i in range(1, k + 1):
            denominator = 2 * n - 2 * i + 1
            sum_terms_list.append(f"1/{denominator}")
        sum_str = " + ".join(sum_terms_list)
        # Add parenthesis for clarity if there are multiple terms
        if k > 1:
            sum_str = f"({sum_str})"
    else:
        # If k=0, no depots are made, and the summation is 0.
        sum_str = "0"

    # Step 3: Print the final equation with all numbers plugged in.
    # The general formula is: (n-k)*100 - (m - 100 * Sum) * (2*(n-k)-1)
    print("The final answer is expressed by the equation:")
    final_equation = f"({n}-{k})*100 - ({m} - 100 * {sum_str}) * (2*({n}-{k})-1)"
    print(final_equation)

    # Calculate the final numerical answer for verification
    sum_val = 0.0
    for i in range(1, k + 1):
        sum_val += 1.0 / (2 * n - 2 * i + 1)
    
    water_left = (n - k) * 100 - (m - 100 * sum_val) * (2 * (n - k) - 1)
    # The final answer is the value of 'water_left'
    
# --- User Input ---
# You can change these values to test different scenarios
n_barrels = 3 # Represents 300 liters
distance_m = 80 # kilometers

solve_horse_problem(n_barrels, distance_m)

# Example 2: n=3, m=15
# solve_horse_problem(3, 15)
# Expected Output: (3-0)*100 - (15 - 100 * 0) * (2*(3-0)-1)

# Example 3: n=2, m=40
# solve_horse_problem(2, 40)
# Expected Output: (2-1)*100 - (40 - 100 * 1/3) * (2*(2-1)-1)