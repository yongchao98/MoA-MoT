import sympy

def solve_transience_condition():
    """
    This function calculates the infimum of c for which the given Markov chain is transient.
    It uses a diffusion approximation method.
    """
    # Define the symbolic variables for state k and parameter c
    k = sympy.Symbol('k')
    c = sympy.Symbol('c')

    # Define the transition probabilities for large k as given in the problem
    P_k_minus_2 = sympy.Rational(1, 4)
    P_k_plus_2 = sympy.Rational(1, 4)
    P_k_minus_1 = sympy.Rational(1, 4) - c/k
    P_k_plus_1 = sympy.Rational(1, 4) + c/k

    # The possible jumps from state k are -2, -1, 1, 2
    jumps = [-2, -1, 1, 2]
    # The corresponding probabilities
    probs = [P_k_minus_2, P_k_minus_1, P_k_plus_1, P_k_plus_2]

    # Step 1: Calculate the mean jump (drift) from state k, mu_k
    mu_k = sum(jump * prob for jump, prob in zip(jumps, probs))
    mu_k = sympy.simplify(mu_k)

    # For large k, the drift mu_k behaves like A/k. We find A by taking the limit.
    A = sympy.limit(k * mu_k, k, sympy.oo)

    # Step 2: Calculate the expected squared jump
    E_jump_sq = sum(jump**2 * prob for jump, prob in zip(jumps, probs))
    E_jump_sq = sympy.simplify(E_jump_sq)

    # For large k, mu_k^2 goes to 0, so the variance sigma^2 is the limit of the expected squared jump.
    sigma_sq = sympy.limit(E_jump_sq, k, sympy.oo)
    
    print("Analysis using Diffusion Approximation:")
    print("-" * 35)
    print(f"1. The mean drift for large k is of the form A/k.")
    print(f"   Calculated mean drift mu_k = {mu_k}")
    print(f"   The drift coefficient A = lim_{{k->oo}} (k * mu_k) = {A}\n")
    
    print(f"2. The variance of the jump for large k is a constant sigma^2.")
    print(f"   Calculated expected squared jump E[(jump)^2] = {E_jump_sq}")
    print(f"   The limiting variance sigma^2 = lim_{{k->oo}} E[(jump)^2] = {sigma_sq}\n")
    
    # Step 3 & 4: Apply the transience condition A > sigma^2 / 2 and solve for c.
    # The chain is transient if A > sigma^2 / 2.
    # The critical value of c is found when A = sigma^2 / 2.
    critical_c_eq = sympy.Eq(A, sigma_sq / 2)
    
    # Solve the equation for c to find the infimum value.
    infimum_c = sympy.solve(critical_c_eq, c)[0]

    print(f"3. The condition for the Markov chain to be transient is A > sigma^2 / 2.")
    # The equation being solved is A = sigma^2 / 2
    A_coeff = A.coeff(c)
    rhs_val = sigma_sq / 2
    
    print(f"   Substituting A and sigma^2, the boundary equation is:")
    print(f"   {A_coeff} * c = {rhs_val}")
    print(f"   This means transience occurs when c > {infimum_c}.\n")
    
    print(f"4. The infimum of the set of c for which the chain is transient is the boundary value.")
    print(f"   Infimum c = {infimum_c}")

solve_transience_condition()