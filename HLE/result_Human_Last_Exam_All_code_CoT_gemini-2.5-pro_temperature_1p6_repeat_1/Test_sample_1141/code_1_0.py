import math

def solve_cutting_problem():
    """
    Calculates the maximum value obtainable by cutting a steel billet
    according to the user-provided (flawed) optimization problem formulation.

    The analysis shows that due to the non-overlapping constraints, the optimal
    strategy is to exclusively cut T1 cubes.
    """

    print("Analyzing the optimal cutting strategy based on the provided rules.")
    print("-----------------------------------------------------------------\n")

    # Strategy 1: Prioritize B2 Balls (Value 150)
    # B2-B2 constraint: squared distance >= 64, or distance >= 8.
    # B2 center ranges: 4<=x<=28, 4<=y<=18, z=4.
    num_b2_x = (28 - 4) // 8 + 1
    num_b2_y = (18 - 4) // 8 + 1
    total_b2 = num_b2_x * num_b2_y
    value_b2_strategy = total_b2 * 150

    print(f"Strategy A (B2 balls first):")
    print(f"The rules allow for {total_b2} B2 balls to be placed.")
    print(f"This yields a value of {total_b2} * 150 = {value_b2_strategy}.")
    print("However, the T1-B2 non-overlapping rule (|z_t1 - z_b2| >= 5) makes it impossible")
    print("to place any T1 cubes if B2 balls are present, as |z_t1 - 4| can be at most 3.")
    print("The remaining value would come from B1 balls (price 1), which is not optimal.\n")

    # Strategy 2: Prioritize T1 Cubes (Value 5)
    # T1-T1 constraint: min(|dx|,|dy|,|dz|) >= 2.
    # This means centers must be spaced by at least 2 units on each axis.
    # T1 center ranges: 1<=x<=31, 1<=y<=21, 1<=z<=7.
    
    # Calculate number of T1s that can be placed on a grid with step 2
    x_coords = range(1, 31 + 1, 2)
    y_coords = range(1, 21 + 1, 2)
    z_coords = range(1, 7 + 1, 2)

    num_t1_x = len(x_coords)
    num_t1_y = len(y_coords)
    num_t1_z = len(z_coords)

    total_t1 = num_t1_x * num_t1_y * num_t1_z
    t1_price = 5
    value_t1_strategy = total_t1 * t1_price

    print(f"Strategy B (T1 cubes first):")
    print("This strategy forbids any B2 balls.")
    print("T1 cubes are placed on a regular grid based on their non-overlapping rule.")
    print(f"Number of T1 cubes along x-axis: {num_t1_x}")
    print(f"Number of T1 cubes along y-axis: {num_t1_y}")
    print(f"Number of T1 cubes along z-axis: {num_t1_z}")
    print(f"Total T1 cubes: {num_t1_x} * {num_t1_y} * {num_t1_z} = {total_t1}\n")

    print("Final Calculation:")
    print("------------------")
    print(f"The total value from Strategy B is significantly higher.")
    print(f"The final equation for the maximum value is:")
    print(f"Total Pieces * Price per Piece = Total Value")
    print(f"{total_t1} * {t1_price} = {value_t1_strategy}")
    
    # As the number of pieces is large, we won't print all of them.
    # But the solution consists of placing a T1 cube at each of these grid points.
    print("\nThe solution is to place a T1 cube at each center (x, y, z) where:")
    print("x is in {1, 3, 5, ..., 31}")
    print("y is in {1, 3, 5, ..., 21}")
    print("z is in {1, 3, 5, 7}")

    return value_t1_strategy

# --- Main Execution ---
if __name__ == "__main__":
    max_value = solve_cutting_problem()
    # The final answer in the requested format
    print(f"\n<<<Answering the questions directly:>>>")
    print(f"Is my problem formulation correct? No, the T1 non-overlapping constraints are flawed.")
    print(f"The highest valid solution based on your formulation is: {max_value}")
    # Final return for the platform
    # print(f"<<<{max_value}>>>")