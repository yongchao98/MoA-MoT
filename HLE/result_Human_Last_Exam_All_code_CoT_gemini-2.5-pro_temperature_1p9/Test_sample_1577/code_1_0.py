import math

def calculate_degeneracy(n, m):
    """
    Calculates the ground state degeneracy of the toric code
    with n smooth and m rough holes based on the formula 2^(n+m-2).

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.
    """
    if n <= 0 or m <= 0:
        print("Note: The formula 2^(n+m-2) is typically derived for cases where n > 0 and m > 0.")
        # Other formulas apply for cases with only one boundary type.
        # For n>0, m=0: GSD = 2^(n-1)
        # For m>0, n=0: GSD = 2^(m-1)
        # The provided multiple-choice options seem to rely on the simplified n,m > 0 case.
        if n > 0 and m == 0:
            exponent = n-1
            print(f"For n={n}, m={m}, a more specific formula gives 2^({n}-1) = {2**exponent}.")
            return
        elif m > 0 and n == 0:
            exponent = m-1
            print(f"For n={n}, m={m}, a more specific formula gives 2^({m}-1) = {2**exponent}.")
            return
        else: # n=0, m=0 (sphere)
            print("For n=0, m=0 (a sphere), the degeneracy is 1.")
            return

    # Calculate based on the formula matching the multiple-choice option
    exponent = n + m - 2
    result = 2**exponent

    print("The ground space degeneracy is given by the formula: 2^(n + m - 2)")
    print(f"For n = {n} and m = {m}:")
    # Output the final equation with the numbers substituted.
    print(f"2^({n} + {m} - 2) = 2^{exponent} = {result}")

# Example usage with some values for n and m.
n_smooth = 4
m_rough = 3
calculate_degeneracy(n_smooth, m_rough)