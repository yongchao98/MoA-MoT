def solve_cutting_optimization():
    """
    This function analyzes the provided cutting stock problem and calculates the maximum possible value based on the given constraints.

    The plan is as follows:
    1.  Analyze the problem constraints. It's identified that the non-overlapping rule between B2 balls and T1 cubes is so strict that it prevents them from co-existing in the billet. This simplifies the problem to choosing between a cutting plan with B2s or one with T1s.
    2.  Calculate the maximum value for a T1-based strategy. Since T1s (value 5) are more valuable than B1s (value 1) and their packing constraints relative to each other are similar, filling the billet with T1s is the optimal approach when B2s are excluded.
    3.  Calculate the maximum value for a B2-based strategy.
    4.  Compare the results from the two strategies to find the global maximum value.
    """

    print("Step 1: Analyzing the problem and constraints.")
    print("Based on the provided rules, a B2 piece cannot be placed in the same billet as a T1 piece.")
    print("This is because a B2's z-center must be 4, and the T1-to-B2 constraint requires a z-distance of at least 5, which is impossible within the billet's height.")
    print("Therefore, we must compare two separate scenarios: a billet with B2s and B1s, or a billet with T1s and B1s.\n")

    # --- Scenario 1: Maximize T1 cubes (value 5) and B1 balls (value 1) ---
    print("Step 2: Calculating value for a T1-based cutting plan.")
    print("Because T1 cubes have a much higher value (5) than B1 balls (1), we prioritize filling the billet with T1s.")
    
    # The constraint for placing two T1s is that their centers must be separated by at least 2 units (1cm) in ALL three dimensions (x, y, and z).
    # We find the number of positions on the grid that satisfy this.
    # T1 center valid ranges: x in [1, 31], y in [1, 21], z in [1, 7]
    num_x_positions_t1 = len(range(1, 31 + 1, 2))
    num_y_positions_t1 = len(range(1, 21 + 1, 2))
    num_z_positions_t1 = len(range(1, 7 + 1, 2))
    
    num_t1 = num_x_positions_t1 * num_y_positions_t1 * num_z_positions_t1
    price_t1 = 5
    value_scenario1 = num_t1 * price_t1
    print(f"Under the given constraints, we can fit {num_x_positions_t1} pieces along the x-axis, {num_y_positions_t1} along the y-axis, and {num_z_positions_t1} along the z-axis.")
    print(f"Total T1 pieces = {num_t1}")
    print(f"Maximum value from this scenario = {num_t1} * {price_t1} = {value_scenario1}\n")

    # --- Scenario 2: Maximize B2 balls (value 150) ---
    print("Step 3: Calculating value for a B2-based cutting plan.")
    
    # The constraint for placing two B2s is that their centers must be at least 8 units (4cm) apart.
    # B2 center valid ranges: x in [4, 28], y in [4, 18], z = 4
    # We place B2s on a grid with a spacing of 8 to meet the distance requirement.
    num_x_positions_b2 = len(range(4, 28 + 1, 8))
    num_y_positions_b2 = len(range(4, 18 + 1, 8))
    
    num_b2 = num_x_positions_b2 * num_y_positions_b2
    price_b2 = 150
    value_scenario2 = num_b2 * price_b2
    print(f"We can fit {num_x_positions_b2} pieces along the x-axis and {num_y_positions_b2} along the y-axis.")
    print(f"Total B2 pieces = {num_b2}")
    print(f"Value from B2 pieces alone = {num_b2} * {price_b2} = {value_scenario2}")
    print("Adding B1 balls of value 1 each to the remaining space will not be enough to surpass the value from Scenario 1.\n")
    
    # --- Conclusion ---
    print("Step 4: Comparing scenarios and concluding.")
    if value_scenario1 > value_scenario2:
        print(f"Comparing the two scenarios ({value_scenario1} vs {value_scenario2}), the T1-only plan is optimal.")
        print("\nThe final equation for the maximum value is:")
        print(f"{num_t1} * {price_t1} = {value_scenario1}")
        return value_scenario1
    else:
        print(f"Comparing the two scenarios ({value_scenario1} vs {value_scenario2}), the B2-based plan is optimal.")
        print("\nThe final equation for the maximum value is:")
        print(f"{num_b2} * {price_b2} = {value_scenario2}")
        return value_scenario2


if __name__ == '__main__':
    solve_cutting_optimization()