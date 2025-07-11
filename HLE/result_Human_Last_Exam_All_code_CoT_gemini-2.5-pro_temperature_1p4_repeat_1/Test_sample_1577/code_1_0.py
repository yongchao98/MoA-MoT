def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    on a plane with n smooth holes and m rough holes.

    The formula is GSD = 2^(n + m - 2), assuming n>=1 and m>=1.
    """
    # The problem asks to choose a formula, which is 2^(n + m - 2).
    # We will demonstrate this formula with example values for n and m.
    # We assume n and m are at least 1, as having zero holes of one type
    # would change the formula, and this form matches one of the options.

    if n <= 0 or m <= 0:
        print("This formula assumes n > 0 and m > 0.")
        # For n=1, m>0, the GSD is 2^(m-1).
        # For m=1, n>0, the GSD is 2^(n-1).
        # The general formula is 2^(max(0, n-1) + max(0, m-1)).
        # However, 2^(n+m-2) is one of the choices, suggesting n>=1 and m>=1 is assumed.
        # We will proceed with the calculation to demonstrate the formula from the answer choice.

    exponent = n + m - 2
    result = 2**exponent

    print(f"The ground space degeneracy (GSD) is calculated using the formula: 2^(n + m - 2)")
    print(f"For n = {n} (smooth holes) and m = {m} (rough holes):")
    # The final output prints the equation with all the numbers.
    print(f"GSD = 2^({n} + {m} - 2) = 2^{exponent} = {result}")

# Example usage with n=3 smooth holes and m=2 rough holes.
calculate_toric_code_gsd(n=3, m=2)