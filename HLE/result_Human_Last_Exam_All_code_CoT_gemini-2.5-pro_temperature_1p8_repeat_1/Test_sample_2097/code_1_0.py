import math
from fractions import Fraction

def laguerre_L_at_1(poly_n, alpha):
    """
    Calculates the value of the associated Laguerre polynomial L_n^(alpha)(x) at x=1.
    The definition is L_n^(alpha)(x) = sum_{i=0 to n} C(n+alpha, n-i) * (-x)^i / i!
    """
    total = Fraction(0)
    for i in range(poly_n + 1):
        # Using fractions for precision
        try:
            term = Fraction(math.comb(poly_n + alpha, poly_n - i) * ((-1)**i), math.factorial(i))
            total += term
        except ValueError:
            # Handles cases where comb fails for negative inputs, etc.
            # which shouldn't happen with the parameters from this problem.
            return float('nan')
    return total

def calculate_Mz1(n):
    """
    Calculates M_z(1) for a given number of spins n.
    """
    # Parameters for the Laguerre polynomial L_{n+1}^{(3n-1)}(1)
    poly_n = n + 1
    alpha = 3 * n - 1

    # Check for validity of alpha, must be > -1
    if alpha <= -1:
        return float('nan')

    L_val = laguerre_L_at_1(poly_n, alpha)

    # Calculate the pre-factor
    # Mz = (-2/pi)^n * (n+1)/n^n * L_val
    # We will calculate the rational part first, then the pi part
    
    # Numerator of the rational part
    num = ((-2)**n) * (n + 1) * L_val.numerator
    # Denominator of the rational part
    den = (n**n) * L_val.denominator
    
    # Keep the final value as a fraction to display with pi
    final_fraction = Fraction(num, den)

    # Numerical value for comparison
    numerical_val = ((-2 / math.pi)**n) * ((n + 1) / (n**n)) * float(L_val)

    return final_fraction, numerical_val

def main():
    """
    Main function to find the minimum M_z(1).
    """
    min_Mz = float('inf')
    n_min = -1
    min_fraction = Fraction(0)
    results = []

    print("Calculating M_z(1) for different values of n:")
    for n in range(1, 11):
        fraction_val, numerical_val = calculate_Mz1(n)
        results.append((n, fraction_val, numerical_val))
        print(f"n = {n}: M_z(1) = {numerical_val:.4f}")
        if numerical_val < min_Mz:
            min_Mz = numerical_val
            n_min = n
            min_fraction = fraction_val
    
    print("\n------------------------------------------------------")
    print(f"The minimum magnetization M_z(1) occurs at n = {n_min}.")
    print(f"The minimum value is {min_Mz:.7f}")
    
    num = min_fraction.numerator
    den = min_fraction.denominator
    pi_power = n_min

    # Simplify the fraction part of the result for nice printing
    # Example for n=3: num = -29380, den = 81, pi_power = 3
    # which is -29380 / (81 * pi^3)

    print("\nThe exact expression for the minimum magnetization is:")
    print(f"M_z(1, n={n_min}) = {num} / ({den} * pi^{pi_power})")
    print("\nBreaking down the calculation for the final answer:")

    n = n_min
    poly_n = n + 1
    alpha = 3 * n - 1
    L_val_frac = laguerre_L_at_1(poly_n, alpha)
    prefactor_num = ((-2)**n) * (n + 1)
    prefactor_den_n_part = n**n
    
    print(f"n_min = {n}")
    print(f"M_z(1, {n}) = (-2/pi)^{n} * ({n}+1)/{n}^{n} * L_{{{n}+1}}^{({3*n}-1)}(1)")
    print(f"L_{poly_n}^{({alpha})}(1) = {L_val_frac.numerator}/{L_val_frac.denominator}")
    print(f"M_z(1, {n}) = ({prefactor_num}/pi^{n}) * (1/{prefactor_den_n_part}) * ({L_val_frac.numerator}/{L_val_frac.denominator})")
    
    final_num = prefactor_num * L_val_frac.numerator
    final_den = prefactor_den_n_part * L_val_frac.denominator
    
    print(f"M_z(1, {n}) = {final_num} / ({final_den} * pi^{n})")

    # To show the equation for the final answer more clearly, print each number
    print("\nThe final equation is:")
    print(f"M_z(1) = {final_num} / ({final_den} * \u03c0^{n})")

if __name__ == "__main__":
    main()