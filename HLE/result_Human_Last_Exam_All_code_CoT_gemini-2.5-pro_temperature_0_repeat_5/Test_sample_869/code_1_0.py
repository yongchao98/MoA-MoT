import math

def calculate_probability():
    """
    Calculates and prints the total, favorable, and probability of the item distribution.
    """
    # Define constants from the problem
    N = 25  # Total items
    T = 5   # Number of types / individuals
    k = 5   # Items per type / per individual

    # Step 1: Calculate the total number of ways to distribute the items (S)
    # S = 25! / (5!)^5
    try:
        s_numerator = math.factorial(N)
        s_denominator = math.factorial(k) ** T
        S = s_numerator // s_denominator
    except OverflowError:
        print("Error: Calculation of S resulted in a number too large to handle.")
        return

    # Step 2: Calculate the number of favorable distributions (F)
    # F = 5! * (C(5, 3))^5 * (10! / 2!^5)
    try:
        # Number of ways to assign specialty types
        f_assignments = math.factorial(T)

        # C(5, 3) = 5! / (3! * 2!)
        comb_5_3 = math.factorial(k) // (math.factorial(3) * math.factorial(k - 3))

        # Number of ways to place the core specialty items
        f_placements = comb_5_3 ** T

        # Number of ways to arrange the remaining 10 items (2 of each of the 5 types)
        remaining_items = N - T * 3
        items_per_type_remaining = k - 3
        f_arrangements = math.factorial(remaining_items) // (math.factorial(items_per_type_remaining) ** T)

        F = f_assignments * f_placements * f_arrangements
    except OverflowError:
        print("Error: Calculation of F resulted in a number too large to handle.")
        return

    # Step 3: Calculate the probability P = F / S
    P = F / S

    # Output the results as per the user's request
    print(f"The total number of ways to distribute the items is:")
    print(f"S = {S}")
    print("\nThe number of favorable distributions is:")
    print(f"F = {F}")
    print("\nThe probability is P = F / S:")
    print(f"P = {F} / {S}")
    print(f"P ≈ {P}")
    
    # For verification, the exact fraction is 1250000 / 572605033
    exact_numerator = 1250000
    exact_denominator = 572605033
    print(f"\nThe exact probability as a fraction is: {exact_numerator}/{exact_denominator}")


if __name__ == "__main__":
    calculate_probability()