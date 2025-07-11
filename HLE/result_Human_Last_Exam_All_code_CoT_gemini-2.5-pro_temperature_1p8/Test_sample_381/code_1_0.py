import math

def get_upper_bound_factor():
    """
    This function derives and prints the upper-bound for ||B * Q_{0,M}||_oo
    as a factor of sqrt(N), based on the provided mathematical context.
    """
    
    # From the mathematical background provided, especially results concerning the
    # ergodicity of products of substochastic matrices (like D*P), we can find a bound
    # for the spectral norm of B * Q_{0,M}. A relevant lemma, which can be found in the
    # literature referenced by the text, states:
    # ||B * Q_{0,M}||_2 <= [product_{t=0 to M} (1 - c * delta_t)]^(1/2)

    # The problem defines beta_M = product_{t=0 to M} (1 - c * delta_t).
    # Using this definition, the inequality for the spectral norm becomes:
    # ||B * Q_{0,M}||_2 <= (beta_M)^(1/2)

    # The goal is to find a bound for the infinity norm, ||.||_oo.
    # We use the standard matrix norm inequality that for a matrix A with N columns,
    # ||A||_oo <= sqrt(N) * ||A||_2.

    # Applying this inequality to our matrix A = B * Q_{0,M} (which has N columns), we get:
    # ||B * Q_{0,M}||_oo <= sqrt(N) * ||B * Q_{0,M}||_2

    # Combining the two inequalities, we arrive at the final bound:
    # ||B * Q_{0,M}||_oo <= sqrt(N) * (beta_M)^(1/2)

    # The question asks for the upper bound expressed as a factor of sqrt(N).
    # From the expression above, this factor is (beta_M)^(1/2).
    
    # We will now print this factor. The formula contains the numbers 1 and 2.
    factor_formula = "(beta_M)^(1/2)"
    
    print(f"The upper-bound for ||B * Q_{0, M}||_oo can be expressed as F * sqrt(N).")
    print(f"The factor F is derived to be: {factor_formula}")
    print("This can also be written as sqrt(beta_M).")
    print("\nIn this final equation for the factor, the numbers are:")
    print("beta_M raised to the power of (1 / 2)")

get_upper_bound_factor()