import math
from fractions import Fraction

def calculate_cn_prefactor(n):
    """
    Calculates the system-independent prefactor c_n for the fully 
    f-connected Ree-Hoover diagram in the virial expansion.

    The formula for c_n is derived from the Ree-Hoover representation of the
    virial coefficients:
    c_n = - (n-1) / n! * (-1)^v(K_n)
    where v(K_n) is the cyclomatic number of the complete graph K_n,
    given by v(K_n) = (n-1)(n-2)/2.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    # Print the general formula
    print(f"The formula for the prefactor c_n is:")
    print(f"c_n = - (n-1) / n! * (-1)^((n-1)*(n-2)/2)")
    print("-" * 30)
    
    # Substitute the given value of n
    print(f"For n = {n}, the calculation is:")
    print(f"c_{n} = - ({n}-1) / {n}! * (-1)^(({n}-1)*({n}-2)/2)")

    # Calculate the components of the formula
    term_n_minus_1 = n - 1
    term_n_factorial = math.factorial(n)
    
    exp_term1 = n - 1
    exp_term2 = n - 2
    exponent = (exp_term1 * exp_term2) // 2

    print(f"c_{n} = - {term_n_minus_1} / {term_n_factorial} * (-1)^({exp_term1}*{exp_term2}/2)")
    print(f"c_{n} = - {term_n_minus_1} / {term_n_factorial} * (-1)^{exponent}")

    # Calculate the sign
    sign = (-1)**exponent
    print(f"c_{n} = - {term_n_minus_1} / {term_n_factorial} * ({sign})")
    
    # Use Fractions for a clean representation of the result
    coeff_frac = Fraction(term_n_minus_1, term_n_factorial)
    final_frac = -1 * coeff_frac * sign
    
    print(f"c_{n} = -({coeff_frac}) * ({sign})")
    print(f"c_{n} = {final_frac}")


# We can demonstrate the calculation for a specific n, for example n=5.
n_value = 5
calculate_cn_prefactor(n_value)