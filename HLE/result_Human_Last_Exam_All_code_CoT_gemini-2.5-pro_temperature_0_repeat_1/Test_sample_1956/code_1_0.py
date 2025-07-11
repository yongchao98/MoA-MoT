import math

def solve_game_positions():
    """
    Calculates the number of starting positions where the bitwise xor value
    of the piles' Grundy values is equal to one or two.
    """
    # The problem states n > 200 and t > 0.
    # We will use example values n = 201 and t = 10.
    n = 201
    t = 10

    print(f"Solving for n = {n} and t = {t}\n")

    # Step 1: Determine the Grundy values (g-values).
    # Let G(a) be the g-value for a pile of size 2^a.
    # The sequence G(a) for a >= 1 is periodic with period 4: (1, 2, 3, 0, 1, 2, 3, 0, ...).
    # G(a) = 1 if a = 1 (mod 4)
    # G(a) = 2 if a = 2 (mod 4)
    # G(a) = 3 if a = 3 (mod 4)
    # G(a) = 0 if a = 0 (mod 4)

    # Step 2: Count the number of available choices for a_i for each g-value.
    # The range for a_i is 1 <= a_i <= 4t+2.
    
    # N1: g-value is 1 (a = 1, 5, ..., 4t+1)
    N1 = t + 1
    # N2: g-value is 2 (a = 2, 6, ..., 4t+2)
    N2 = t + 1
    # N3: g-value is 3 (a = 3, 7, ..., 4t-1)
    N3 = t
    # N0: g-value is 0 (a = 4, 8, ..., 4t)
    N0 = t
    
    print("Number of choices for a_i for each Grundy value:")
    print(f"N0 (g-value 0): {N0}")
    print(f"N1 (g-value 1): {N1}")
    print(f"N2 (g-value 2): {N2}")
    print(f"N3 (g-value 3): {N3}\n")

    # Step 3: Use the derived combinatorial formula to find the total count.
    # The number of positions S is given by the formula:
    # S = 2^(n-1) * ((2*t + 1)^n - (-1)^n)
    
    # Let's calculate the components of the formula.
    term1_base = 2 * t + 1
    term2_base = -1
    power_of_2 = n - 1

    print("The formula for the number of positions is: 2**(n-1) * ((2*t + 1)**n - (-1)**n)")
    print("Substituting the values of n and t:")
    # Using f-string to display the equation with the numbers
    final_equation = f"2**({power_of_2}) * (({term1_base})**{n} - ({term2_base})**{n})"
    print(f"Equation: {final_equation}\n")

    # Step 4: Calculate the final result.
    # Python's integers handle arbitrary size, so we can compute this directly.
    result = (2**power_of_2) * (term1_base**n - term2_base**n)

    print(f"The total number of starting positions is: {result}")

solve_game_positions()