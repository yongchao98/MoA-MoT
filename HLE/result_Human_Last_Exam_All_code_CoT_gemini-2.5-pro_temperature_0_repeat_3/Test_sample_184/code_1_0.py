import sympy

def solve_sum():
    """
    Calculates the sum of 1/n^2 for all n in S, where S is the set of positive
    integers n such that no prime divides n with multiplicity 1, 2, or 5.

    The derivation shows the sum is equal to zeta(6) * zeta(8).
    """

    # The sum is equivalent to zeta(6) * zeta(8)
    # Get the values for zeta(6) and zeta(8)
    # zeta(6) = pi^6 / 945
    # zeta(8) = pi^8 / 9450
    
    # We can use sympy to verify, but for the final print, we'll use the known fractions.
    # z6 = sympy.zeta(6)
    # z8 = sympy.zeta(8)
    
    den_z6 = 945
    den_z8 = 9450
    
    # The product is (pi^6 / 945) * (pi^8 / 9450)
    final_power_of_pi = 6 + 8
    final_denominator = den_z6 * den_z8
    
    # Print the steps of the final calculation
    print("The sum is derived to be the product of two Riemann zeta function values: zeta(6) and zeta(8).")
    print(f"Sum = zeta(6) * zeta(8)")
    print(f"      = (pi^6 / {den_z6}) * (pi^8 / {den_z8})")
    print(f"      = pi^({6 + 8}) / ({den_z6} * {den_z8})")
    print(f"      = pi^{final_power_of_pi} / {final_denominator}")

solve_sum()
<<<pi^14 / 8930250>>>