import sympy

def solve_sum():
    """
    Calculates the sum of 1/n^2 for all positive integers n where no prime
    divides n with multiplicity 1, 2, or 5.
    The solution is expressed as a rational number times a power of pi.
    """
    print("The problem is to evaluate the sum of 1/n^2 for all n in a set S.")
    print("S is the set of positive integers n such that for any prime p, the exponent of p in the factorization of n is not 1, 2, or 5.")
    print("This sum can be written as an Euler product, which simplifies to the product of two Riemann zeta functions: zeta(6) * zeta(8).\n")

    # We use the formula for zeta(2n) involving Bernoulli numbers:
    # zeta(2n) = (-1)^(n+1) * B_{2n} * (2*pi)^(2n) / (2 * (2n)!)
    
    # Calculate zeta(6)
    n_6 = 3
    B_6 = sympy.bernoulli(2 * n_6)  # B_6 = 1/42
    zeta_6_expr = (-1)**(n_6 + 1) * B_6 * (2 * sympy.pi)**(2 * n_6) / (2 * sympy.factorial(2 * n_6))
    zeta_6_val = sympy.simplify(zeta_6_expr)

    # Calculate zeta(8)
    n_8 = 4
    B_8 = sympy.bernoulli(2 * n_8)  # B_8 = -1/30
    zeta_8_expr = (-1)**(n_8 + 1) * B_8 * (2 * sympy.pi)**(2 * n_8) / (2 * sympy.factorial(2 * n_8))
    zeta_8_val = sympy.simplify(zeta_8_expr)

    # The final sum is the product of zeta(6) and zeta(8)
    total_sum = zeta_6_val * zeta_8_val

    # Extract the components for printing the final equation
    c6 = zeta_6_val / (sympy.pi**6)
    c8 = zeta_8_val / (sympy.pi**8)
    c6_num, c6_den = c6.as_numer_denom()
    c8_num, c8_den = c8.as_numer_denom()
    
    final_coeff = total_sum / (sympy.pi**14)
    fc_num, fc_den = final_coeff.as_numer_denom()

    print("Step 1: Calculate zeta(6)")
    print(f"zeta(6) = {zeta_6_val}\n")
    
    print("Step 2: Calculate zeta(8)")
    print(f"zeta(8) = {zeta_8_val}\n")

    print("Step 3: Multiply the results to get the final sum.")
    print("The final equation is:")
    print(f"Sum = zeta(6) * zeta(8)")
    print(f"    = ({c6_num}/{c6_den} * pi^6) * ({c8_num}/{c8_den} * pi^8)")
    print(f"    = ({c6_num * c8_num}/{c6_den * c8_den}) * pi^(6 + 8)")
    print(f"    = {fc_num}/{fc_den} * pi^14\n")

    print(f"The final answer is the rational number {final_coeff} times pi to the power of 14.")

if __name__ == '__main__':
    solve_sum()