import sympy as sp

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for which the given Markov chain is transient.
    """
    # Define symbols for symbolic computation
    k = sp.Symbol('k')
    c = sp.Symbol('c')

    print("Step 1: Define the transition probabilities for jumps from state k.")
    # Jumps are -2, -1, 1, 2. The probabilities P_k,k+j are denoted P[j].
    P = {
        -2: sp.Rational(1, 4),
        -1: sp.Rational(1, 4) - c/k,
        1:  sp.Rational(1, 4) + c/k,
        2:  sp.Rational(1, 4)
    }
    print("P(jump=-2) =", P[-2])
    print("P(jump=-1) =", P[-1])
    print("P(jump=+1) =", P[1])
    print("P(jump=+2) =", P[2])
    print("-" * 30)

    print("Step 2: Calculate the expected displacement (drift) mu_k.")
    mu_k = sp.sympify(0)
    for jump, prob in P.items():
        mu_k += jump * prob
    mu_k = sp.simplify(mu_k)
    print(f"mu_k = E[Delta X] = {mu_k}")
    print("The drift is of the form a/k. We can identify 'a' by calculating lim_{k->inf} (k * mu_k).")
    a = sp.limit(k * mu_k, k, sp.oo)
    print(f"The drift coefficient 'a' is: {a}")
    print("-" * 30)

    print("Step 3: Calculate the limiting variance of the displacement, sigma^2.")
    E_sq_jump = sp.sympify(0)
    for jump, prob in P.items():
        E_sq_jump += jump**2 * prob
    E_sq_jump = sp.simplify(E_sq_jump)
    print(f"E[(Delta X)^2] = {E_sq_jump}")
    
    variance_k = E_sq_jump - mu_k**2
    sigma_sq = sp.limit(variance_k, k, sp.oo)
    print(f"The limiting variance sigma^2 = lim_{k->inf} Var(Delta X) = {sigma_sq}")
    print("-" * 30)
    
    print("Step 4: Apply the transience criterion.")
    print("A random walk with drift mu_k ~ a/k and limiting variance sigma^2 is transient if a > sigma^2 / 2.")
    
    # The criterion for transience is a > sigma^2 / 2
    # We substitute the found values into this inequality.
    # a = 2*c
    # sigma^2 = 5/2
    lhs_expr = a
    rhs_expr = sigma_sq / 2
    
    print("\nSubstituting the values into the inequality a > sigma^2 / 2:")
    print(f"The term 'a' is: {lhs_expr}")
    print(f"The term 'sigma^2 / 2' is: {rhs_expr}")
    
    print(f"\nThe inequality is: {lhs_expr} > {rhs_expr}")
    
    # Solve the inequality for c
    # In this case, it's 2*c > 5/4, which simplifies to c > 5/8
    infimum_c = rhs_expr / (a/c) # (5/4) / (2c/c) = (5/4) / 2 = 5/8

    print(f"\nSolving for c, we get: c > {infimum_c}")
    print("-" * 30)

    print("Step 5: Determine the infimum.")
    print("The set of c for which the chain is transient is (5/8, infinity).")
    print(f"The infimum of this set is {infimum_c}.")

if __name__ == '__main__':
    solve_markov_chain_transience()
