import math

def solve_radiation_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    (pi/4 half-angle) to the intensity along a line for the given physical setup.
    
    The final expression for the ratio is S = 4 * pi * (4/3 - 7/(6*sqrt(2))).
    """

    # --- Component numbers of the final equation ---
    num_4_pi = 4.0
    num_pi = math.pi
    num_4_over_3 = 4.0 / 3.0
    num_7 = 7.0
    num_6 = 6.0
    num_sqrt_2 = math.sqrt(2.0)
    
    # --- Calculation ---
    term_in_parentheses = num_4_over_3 - num_7 / (num_6 * num_sqrt_2)
    result = num_4_pi * num_pi * term_in_parentheses

    # --- Output the results as requested ---
    print("The final equation for the ratio is S = 4 * pi * (4/3 - 7 / (6 * sqrt(2)))")
    print("\n--- Numerical values in the equation ---")
    print(f"The first number is: {num_4_pi}")
    print(f"The value of pi is: {num_pi}")
    print(f"The fraction 4/3 is: {num_4_over_3}")
    print(f"The number 7 is: {num_7}")
    print(f"The number 6 is: {num_6}")
    print(f"The value of sqrt(2) is: {num_sqrt_2}")

    print("\n--- Final result ---")
    print(f"The maximum achievable ratio is: {result}")

solve_radiation_ratio()