def solve_turbo_snail_problem():
    """
    Calculates the minimum number of attempts Turbo needs to guarantee
    reaching the last row.
    """

    # Step 1: Define the parameters from the problem description.
    num_rows = 2024
    num_cols = 2023
    num_monsters = 2022

    print("--- Problem Analysis ---")
    print(f"The board has {num_cols} columns.")
    print(f"There are {num_monsters} monsters, each in a unique column.")
    print("This means there are 2022 unsafe columns and exactly one safe column.")
    print("")

    # Step 2: Frame the problem as a search.
    # Turbo's goal is to find the single safe column out of all possibilities.
    num_candidates = num_cols
    num_safe_columns = num_cols - num_monsters

    print("--- Strategy for a Guaranteed Win ---")
    print(f"Turbo must find the {num_safe_columns} safe column among {num_candidates} total columns.")
    print("A guaranteed strategy must work for the worst-case monster placement.")
    print("In the worst case, each failed attempt only reveals one unsafe column.")
    print("")

    # Step 3: Calculate the number of attempts needed.
    # To be certain, Turbo must eliminate all other possibilities.
    # The number of columns to eliminate is the number of unsafe columns.
    num_eliminations_needed = num_candidates - num_safe_columns

    # In the worst-case, this requires one failed attempt per column to be eliminated.
    worst_case_failed_attempts = num_eliminations_needed

    # After all unsafe columns are found, one final attempt is needed for the safe column.
    final_successful_attempt = 1

    # The total number of attempts is the sum of the worst-case failures and the final success.
    total_attempts = worst_case_failed_attempts + final_successful_attempt

    print("--- Calculation ---")
    print(f"Number of columns to eliminate = Total Columns - Safe Columns")
    print(f"Number of columns to eliminate = {num_candidates} - {num_safe_columns} = {num_eliminations_needed}")
    print("")
    print(f"Worst-case failed attempts required = {worst_case_failed_attempts}")
    print(f"Final successful attempt required = {final_successful_attempt}")
    print("")
    print("Total guaranteed attempts (n) = (Worst-case failed attempts) + (Final successful attempt)")
    print(f"n = {worst_case_failed_attempts} + {final_successful_attempt}")
    print(f"n = {total_attempts}")

solve_turbo_snail_problem()