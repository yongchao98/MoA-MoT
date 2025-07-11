import sympy as sp

def solve_steady_state_pi0():
    """
    Solves for the steady-state probability pi_0 for the given birth-death process.
    
    The process has rates:
    - Birth rate (i -> i+1): lambda_i = lambda / (i+1)
    - Death rate (i -> i-1): mu_i = mu
    
    The solution is expressed in terms of rho = lambda / mu.
    """
    # Define the necessary symbols for our symbolic computation
    # rho represents lambda / mu
    rho = sp.Symbol('rho', real=True, positive=True)
    # n is the state index
    n = sp.Symbol('n', integer=True, nonnegative=True)

    # Step 1: Find the general expression for pi_n in terms of pi_0.
    # The balance equations are pi_i * lambda_i = pi_{i+1} * mu_{i+1}.
    # pi_{i+1} = pi_i * (lambda_i / mu_{i+1})
    #   lambda_i = lambda / (i+1)
    #   mu_{i+1} = mu
    # pi_{i+1} = pi_i * (lambda / (mu * (i+1)))
    # By induction, this leads to pi_n = pi_0 * (rho**n / n!)
    pi_n_ratio = rho**n / sp.factorial(n)
    
    print("Step 1: Express pi_n in terms of pi_0")
    print(f"The ratio pi_n / pi_0 is given by: {pi_n_ratio}")
    print("-" * 50)

    # Step 2: Use the normalization condition: Sum(pi_n for n=0 to inf) = 1.
    # This implies pi_0 * Sum(pi_n / pi_0 for n=0 to inf) = 1.
    # We calculate the sum of the series.
    series_sum = sp.summation(pi_n_ratio, (n, 0, sp.oo))
    
    print("Step 2: Apply the normalization condition Sum(pi_n) = 1")
    print(f"This requires calculating the sum of pi_n/pi_0 from n=0 to infinity.")
    print(f"Sum({pi_n_ratio}, (n, 0, oo)) = {series_sum}")
    print("-" * 50)
    
    # Step 3: Solve for pi_0.
    # pi_0 * series_sum = 1  =>  pi_0 = 1 / series_sum
    pi_0_solution = 1 / series_sum

    print("Step 3: Solve for pi_0")
    print(f"From pi_0 * ({series_sum}) = 1, we get:")
    print(f"pi_0 = {pi_0_solution}")
    print("\nThe final equation is pi_0 = exp(-1 * rho).")
    print("The numbers in this final equation are:")
    print("The subscript of pi: 0")
    print("The coefficient of rho in the exponent: -1")

# Execute the function to print the derivation and result
solve_steady_state_pi0()