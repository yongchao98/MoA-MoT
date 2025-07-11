import sympy as sp

def solve_problem():
    """
    This function calculates the symbolic expression for phi(n) based on the step-by-step derivation.
    The derivation simplifies the problem under the assumption that the problem intends for M=J/n
    to be on the manifold, which necessitates I1=I2.
    
    The steps are:
    1. Simplify the projection operator based on the tangent space H*1=0.
    2. Use the identity det(exp(A)) = exp(tr(A)).
    3. Calculate the trace of the projected matrix Proj_M(X^-1).
    4. The trace calculation involves finding the inverse of X and summing its elements in a specific way.
       - tr(X^-1) = 2n
       - tr(X^-1 * P) = (4n-2)/n
    5. The trace of the projection is tr(X^-1) - tr(X^-1 * P) = 2n - (4n-2)/n = 2(n-1)^2 / n.
    6. The final result is the exponential of this trace.
    """
    n = sp.Symbol('n')
    
    # Trace of the projected matrix, as derived in the plan.
    trace_H = 2 * (n - 1)**2 / n
    
    # phi(n) is the exponential of the trace.
    phi_n = sp.exp(trace_H)
    
    # Extract the numbers from the final expression to print them as requested.
    # The expression is exp(2 * (n - 1)**2 / n)
    # The numbers are 2, 1, 2, 1
    
    num1 = 2
    num2 = 1
    num3 = 2
    
    print(f"The calculation leads to the following expression for phi(n):")
    # Using sp.pretty_print for a nicer mathematical layout.
    # We construct the equation phi(n) = expression
    equation = sp.Eq(sp.Symbol('phi(n)'), phi_n)
    
    # We want to print the final equation including the numbers, as requested.
    # The format is phi(n) = exp(num1 * (n - num2)**num3 / n)
    
    print(f"phi(n) = exp({num1}*(n - {num2})**{num3} / n)")


solve_problem()