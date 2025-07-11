import sympy

def solve_transience_condition():
    """
    This script finds the infimum of c for which the given Markov chain is transient.
    It uses the criterion for random walks on integers, which relates transience
    to the asymptotic drift and variance of jumps.
    """
    # Define symbols for our calculation
    k, c = sympy.symbols('k c')

    # The jumps from state k are: -2, -1, 1, 2.
    # The corresponding transition probabilities for large k are:
    # P(k -> k-2) = 1/4
    # P(k -> k-1) = 1/4 - c/k
    # P(k -> k+1) = 1/4 + c/k
    # P(k -> k+2) = 1/4
    jumps = [-2, -1, 1, 2]
    probs = [
        sympy.Rational(1, 4),
        sympy.Rational(1, 4) - c/k,
        sympy.Rational(1, 4) + c/k,
        sympy.Rational(1, 4)
    ]

    print("Step 1: Calculate the asymptotic drift mu_k = E[X_{n+1} - k | X_n = k]")
    # The drift is the expected value of the jump size from state k.
    mu_k = sum(jump * prob for jump, prob in zip(jumps, probs))
    mu_k = sympy.simplify(mu_k)

    print(f"The drift mu_k is calculated as: (-2)*({probs[0]}) + (-1)*({probs[1]}) + (1)*({probs[2]}) + (2)*({probs[3]})")
    print(f"Simplified drift: mu_k = {mu_k}")

    # For large k, the drift is of the form A/k. We find A.
    A = sympy.limit(mu_k * k, k, sympy.oo)
    print(f"The drift is of the form A/k, where the coefficient A = {A}.")
    print("-" * 40)

    print("Step 2: Calculate the asymptotic variance of jumps sigma_k^2 = E[(X_{n+1} - k)^2 | X_n = k]")
    # The variance of jumps is the expected value of the squared jump size.
    sigma_k_sq = sum(jump**2 * prob for jump, prob in zip(jumps, probs))
    sigma_k_sq = sympy.simplify(sigma_k_sq)
    
    print(f"The variance sigma_k^2 is calculated as: (-2)^2*({probs[0]}) + (-1)^2*({probs[1]}) + (1)^2*({probs[2]}) + (2)^2*({probs[3]})")
    print(f"Simplified variance: sigma_k^2 = {sigma_k_sq}")
    
    # For large k, the variance converges to a constant B.
    B = sympy.limit(sigma_k_sq, k, sympy.oo)
    print(f"The variance of jumps converges to a constant B = {B}.")
    print("-" * 40)

    print("Step 3: Apply the transience criterion 2*A > B and solve for c.")
    # The criterion for transience is 2*A > B.
    inequality_lhs_expr = 2 * A
    inequality_rhs_expr = B
    
    print(f"The inequality for transience is: 2 * A > B")
    print(f"Substituting the values of A and B:")
    print(f"2 * ({A}) > {B}")
    
    final_lhs = sympy.simplify(inequality_lhs_expr)
    print("This simplifies to the final inequality:")
    print(f"{final_lhs} > {inequality_rhs_expr}")

    # Solve the inequality for c
    c_coefficient = final_lhs.coeff(c)
    c_infimum = inequality_rhs_expr / c_coefficient
    print(f"\nSolving for c gives: c > {inequality_rhs_expr} / {c_coefficient}")
    print(f"c > {c_infimum}")
    print("-" * 40)
    
    print(f"The condition for the chain to be transient is c > {c_infimum}.")
    print(f"The infimum of the set of such c is the boundary value.")
    final_answer = c_infimum
    print(f"\nFinal Answer: The infimum of c is {final_answer}.")
    
if __name__ == "__main__":
    solve_transience_condition()