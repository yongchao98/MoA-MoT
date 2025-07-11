import math

def solve_bound_problem():
    """
    This function outlines the derivation for the upper-bound and prints the result.
    """
    
    # The plan is to first find a bound for the spectral norm (2-norm)
    # and then convert it to a bound for the infinity norm.

    # 1. The condition `beta > 0` implies that `Sum(delta_t)` over t=0 to infinity is finite.
    #    This means the cumulative perturbation to the chain of stochastic matrices is bounded.

    # 2. Based on matrix perturbation theory, this implies that the spectral norm
    #    of the projected product is bounded by a constant C (assuming a uniform spectral gap
    #    for the class of graphs, making C independent of N).
    #    So, ||B * Q_{0, M}||_2 <= C

    # 3. We use the standard matrix norm inequality that relates the infinity norm and the 2-norm:
    #    ||A||_inf <= sqrt(N) * ||A||_2 for any N x N matrix A.

    # 4. Combining these findings:
    #    ||B * Q_{0, M}||_inf <= sqrt(N) * ||B * Q_{0, M}||_2 <= C * sqrt(N)

    # 5. The question asks for the upper-bound "expressed as a factor of sqrt(N)".
    #    This can be interpreted as the bound on the ratio: ||B * Q_{0, M}||_inf / sqrt(N).
    #    From the previous step, this ratio is bounded by C.
    
    # We will print the final resulting inequality, including the numbers `0`, `M`, and `0.5` (for sqrt).
    
    print("The final inequality for the upper bound is:")
    
    # Define the symbols for printing the equation
    B = "B"
    Q = "Q"
    N = "N"
    C = "C"
    sub_0 = "0"
    sub_M = "M"
    norm_inf = "inf"
    power_half = "0.5" # for sqrt

    # Construct and print the equation string
    # ||B * Q_(0, M)||_inf <= C * N^0.5
    final_equation = f"||{B} * {Q}_({sub_0}, {M})||_{norm_inf} <= {C} * {N}^{power_half}"
    
    print(final_equation)

if __name__ == '__main__':
    solve_bound_problem()