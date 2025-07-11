import numpy as np

def solve_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula.
    
    The problem simplifies under the constraint that M must lie on the manifold,
    which implies I1 = I2. This leads to a simplified tangent space and allows
    for the direct calculation of the trace of the projected matrix.

    The trace of the projection P = Proj_M(X^-1) is found to be:
    Tr(P) = 2(n-1)^2 / n

    The final result phi(n) is exp(Tr(P)).
    """
    
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")
        
    # Step 1: Calculate the trace of the projected matrix P
    # Tr(P) = 2(n-1)^2 / n
    trace_P = 2 * (n - 1)**2 / n
    
    # Step 2: Calculate phi(n) which is exp(Tr(P))
    phi_n = np.exp(trace_P)
    
    # Step 3: Print the results showing each number in the final equation
    print(f"For n = {n}:")
    print(f"The trace of the projected matrix is calculated as:")
    print(f"Tr(P) = 2 * ({n} - 1)^2 / {n} = {trace_P}")
    print(f"The final value is phi({n}) = exp(Tr(P))")
    print(f"phi({n}) = exp({trace_P}) = {phi_n}")
    
# Set n to a value >= 5 as specified in the problem. We choose n=5.
n_value = 5
solve_phi(n_value)
