def solve_turbo_snail_problem():
    """
    This function solves the Turbo the snail problem by calculating the minimum
    number of attempts needed to guarantee reaching the last row.
    """

    # The number of columns on the board.
    num_columns = 2023

    # There are 2022 monsters, and each is in a unique column.
    # This means there are 2022 "unsafe" columns.
    num_unsafe_columns = 2022

    # The core of the problem is for Turbo to find the single "safe" column.
    # A strategy that guarantees success must work for the worst-case scenario.

    # In the worst-case scenario, the monsters are placed such that Turbo must
    # test every single unsafe column before he can identify the safe one.
    # Each time he tests an unsafe column, the attempt fails.
    # An attempt ends when one monster is found, so one failed attempt can, at best,
    # rule out one column.
    max_failed_attempts = num_unsafe_columns

    # After 'max_failed_attempts' have occurred, Turbo has found all 2022 monsters
    # and identified all 2022 unsafe columns. He can now deduce with certainty that
    # the one remaining column is the safe one.

    # His next attempt will be down this known safe column, and it is guaranteed
    # to succeed.
    final_successful_attempt = 1

    # The total number of attempts in this guaranteed strategy is the sum of the
    # maximum possible failures plus the one final successful attempt.
    n = max_failed_attempts + final_successful_attempt

    # The final equation representing the calculation for n.
    print(f"{max_failed_attempts} + {final_successful_attempt} = {n}")

solve_turbo_snail_problem()