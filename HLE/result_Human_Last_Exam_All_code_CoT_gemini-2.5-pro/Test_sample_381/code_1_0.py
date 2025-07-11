import sympy

def solve_task():
    """
    This function calculates and prints the symbolic factor for the upper bound.
    """
    # Define the symbolic variables from the problem description
    # beta: The limit of the product product_{t=0 to k}(1 - c*delta_t) as k -> infinity
    # beta_M: The partial product up to M, i.e., product_{t=0 to M}(1 - c*delta_t)
    # eta_prime: A constant related to the convergence of delta_t, defined as (1-eta)^(b/c)
    beta, beta_M, eta_prime = sympy.symbols("beta beta_M eta'")

    # The upper-bound for the norm ||B * Q_{0,M}||_inf is given in the relevant theorem as:
    # sqrt(N) * (1 - eta')**(-1/2) * beta_M**(-1) * beta
    # We need to find the factor that multiplies sqrt(N).
    
    # Construct the factor expression.
    # It can be written as beta / (beta_M * sqrt(1 - eta_prime))
    factor = beta * (beta_M**-1) * ((1 - eta_prime)**(-sympy.S(1)/2))

    # Print the final result in a clear, readable format.
    # The numbers in the final equation are 1 and -1/2 (or 1/2 for the sqrt).
    # sympy.pprint will display the mathematical formula.
    print("The upper-bound for ||B * Q_{0,M}||_inf is a product of sqrt(N) and a factor.")
    print("The factor is:")
    sympy.pprint(factor)

if __name__ == "__main__":
    solve_task()