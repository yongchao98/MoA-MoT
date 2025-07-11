import math

def solve():
    """
    Calculates the sum of 1/n^2 for all positive integers n such that no prime 
    divides n with multiplicity 1, 2, or 5.
    
    The sum is determined to be zeta(6) * zeta(8).
    zeta(6) = pi^6 / 945
    zeta(8) = pi^8 / 9450
    
    The final sum is (pi^6 / 945) * (pi^8 / 9450) = pi^14 / (945 * 9450).
    """
    
    # The numerators of the rational parts of zeta(6) and zeta(8) are both 1.
    num_zeta6 = 1
    den_zeta6 = 945
    
    num_zeta8 = 1
    den_zeta8 = 9450
    
    # The resulting power of pi is 6 + 8 = 14.
    pi_power = 14
    
    # The final numerator of the rational part is the product of the numerators.
    final_num = num_zeta6 * num_zeta8
    
    # The final denominator of the rational part is the product of the denominators.
    final_den = den_zeta6 * den_zeta8
    
    print("The sum is given by the product of zeta(6) and zeta(8).")
    print(f"zeta(6) = {num_zeta6} * pi**6 / {den_zeta6}")
    print(f"zeta(8) = {num_zeta8} * pi**8 / {den_zeta8}")
    print("The product is (pi**6 / 945) * (pi**8 / 9450)")
    print(f"Resulting in pi**{pi_power} / ({den_zeta6} * {den_zeta8})")
    print("\nFinal Answer:")
    print(f"{final_num} * pi**{pi_power} / {final_den}")

solve()
<<<1 * pi**14 / 8930250>>>