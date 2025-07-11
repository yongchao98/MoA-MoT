import sympy as sp
import math

def solve():
    """
    This function derives the constant 'b' in the asymptotic formula for the
    expected cover and return time C(n) on a uniform random tree.
    The given formula is C(n) ~ b * n^(3/2).
    """

    # Define symbols for our calculation
    n = sp.Symbol('n')
    b = sp.Symbol('b')

    # Step 1: Start with the known asymptotic formula for D_n, the expected
    # distance between two fixed vertices in a uniform random tree.
    # From the literature (e.g., Meir & Moon, 1973), we have: D_n ~ sqrt(pi*n/2)
    D_n_asymptotic = sp.sqrt(sp.pi * n / 2)

    # Step 2: The expected Wiener index E[W(T)] is the sum of expected distances
    # over all pairs of vertices, which is C(n,2) * D_n.
    # For large n, C(n,2) = n*(n-1)/2 ~ n^2/2.
    # E[W(T)] ~ (n^2/2) * D_n
    W_T_asymptotic = (n**2 / 2) * D_n_asymptotic

    # Step 3: The expected cover-and-return time C(n), averaged over starting
    # vertices and trees, is related to the expected Wiener index.
    # C(n) ~ (4/n) * E[W(T)] for large n.
    C_n_asymptotic = (4 / n) * W_T_asymptotic

    # Step 4: Simplify the expression for C(n) to find its asymptotic form.
    C_n_simplified = sp.simplify(C_n_asymptotic)

    # Step 5: Extract the constant 'b' from the resulting expression.
    # The asymptotic expression is of the form b * n^(3/2). The constant 'b' is
    # the coefficient of n^(3/2).
    b_val = C_n_simplified / (n**sp.Rational(3, 2))

    # --- Output the results ---
    print("The derivation for the constant 'b' proceeds as follows:")
    print(f"1. The expected distance between two vertices is D_n ~ {D_n_asymptotic}.")
    print(f"2. The expected Wiener index E[W(T)] ~ (n^2/2) * D_n = {sp.simplify(W_T_asymptotic)}.")
    print(f"3. The expected cover-and-return time C(n) ~ (4/n) * E[W(T)] = {C_n_simplified}.")

    print("\nWe are given the asymptotic relation: C(n) ~ b * n^(3/2).")
    print("By comparing this with our derived formula, we obtain the final equation:")
    
    # Format the output to clearly show all the numbers involved.
    formatted_b = f"sqrt(2)*sqrt(pi)"
    print(f"{formatted_b} * n^(3/2) = b * n^(3/2)")
    
    print(f"\nFrom this equation, we can see that the exact value for 'b' is:")
    print(f"b = sqrt(2*pi)")
    print(f"Numerically, b is approximately {math.sqrt(2 * math.pi):.8f}")

solve()
<<<sqrt(2*pi)>>>