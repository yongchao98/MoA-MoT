import sympy

def solve_bargaining_problem():
    """
    Solves the described bargaining problem by:
    1. Calculating the SPE payoff for Player 1.
    2. Calculating the optimal payoff for Player 1 against the specific Player 2.
    3. Finding the difference and the sum of absolute coefficients of the resulting rational function.
    """
    d = sympy.Symbol('d')

    # Step 1: Calculate Player 1's payoff in the Subgame Perfect Equilibrium (SPE)
    # Let v1 be P1's payoff when P1 proposes.
    # P1 offers x. P2 accepts if 1-x >= d*w2, where w2 is P2's payoff when P2 proposes.
    # P1 offers x such that 1-x = d*w2, so v1 = x = 1 - d*w2.
    # When P2 proposes y, P1 accepts if y >= d*v1. P2 offers y = d*v1.
    # P2's payoff is w2 = 1-y = 1-d*v1.
    # Solving v1 = 1 - d*(1-d*v1) gives v1 = (1-d)/(1-d**2) = 1/(1+d).
    V_SPE = 1 / (1 + d)
    print("--- Step 1: Standard Subgame Perfect Equilibrium ---")
    print(f"Player 1's payoff in the SPE is V_SPE = {V_SPE}\n")

    # Step 2 & 3: Model specific P2 and find P1's optimal payoff
    # P2's decision rule: Accept if payoff(A) > E[payoff(R)].
    # Payoff(A) = 1-x.
    # E[payoff(R)] is P2's expected payoff from the next round, assuming random play.
    # Let W2 be P2's expected payoff when P2 proposes (random play).
    # Let V2 be P2's expected payoff when P1 proposes (random play).
    #
    # When P2 proposes (state W): P2 offers y~U[0,1]. P1 accepts with prob 1/2.
    # W2 = 0.5 * E[1-y] + 0.5 * d * V2. E[1-y]=0.5. So, W2 = 0.25 + 0.5*d*V2.
    # When P1 proposes (state V): P1 offers z~U[0,1] for P2. P2 accepts with prob 1/2.
    # V2 = 0.5 * E[z] + 0.5 * d * W2. E[z]=0.5. So, V2 = 0.25 + 0.5*d*W2.
    
    W2 = sympy.Symbol('W2')
    V2 = sympy.Symbol('V2')
    
    eq1 = sympy.Eq(W2, sympy.Rational(1, 4) + d/2 * V2)
    eq2 = sympy.Eq(V2, sympy.Rational(1, 4) + d/2 * W2)
    
    # Solve the system for W2, which is P2's E[payoff(R)]
    solution = sympy.solve([eq1, eq2], (W2, V2))
    E_R_P2 = solution[W2]

    print("--- Step 2 & 3: Payoff vs. Specific Opponent ---")
    print(f"Player 2's expected payoff from rejecting (assuming random play) is: {sympy.simplify(E_R_P2)}")

    # P2 accepts P1's offer x if 1-x > E_R_P2.
    # P1 wants to maximize x, so P1 sets 1-x = E_R_P2.
    # P1's optimal payoff is x = 1 - E_R_P2.
    V_P1_optimal = 1 - E_R_P2
    V_P1_optimal_simplified = sympy.simplify(V_P1_optimal)
    print(f"Player 1's optimal payoff against this opponent is V_P1_optimal = 1 - ({E_R_P2}) = {V_P1_optimal_simplified}\n")

    # Step 4: Compute the difference
    difference = V_P1_optimal - V_SPE
    
    # Step 5: Express as p(d)/q(d) with integer coefficients
    # sympy.cancel() simplifies the fraction and finds the common denominator
    difference_rational = sympy.cancel(difference)
    
    p_d, q_d = sympy.fraction(difference_rational)

    # To ensure integer coefficients, we can multiply numerator and denominator
    # by the LCM of the denominators of the coefficients.
    p_poly = sympy.Poly(p_d, d)
    q_poly = sympy.Poly(q_d, d)
    lcm_p = sympy.lcm([c.q for c in p_poly.coeffs()])
    lcm_q = sympy.lcm([c.q for c in q_poly.coeffs()])
    common_mult = sympy.lcm(lcm_p, lcm_q)
    
    p_final = sympy.expand(p_d * common_mult)
    q_final = sympy.expand(q_d * common_mult)
    
    # Make leading coefficients positive if possible
    if p_final.as_poly(d).LC() < 0:
        p_final = -p_final
        q_final = -q_final
        
    print("--- Step 4 & 5: Calculate the Difference ---")
    print(f"The difference is V_P1_optimal - V_SPE = ({V_P1_optimal_simplified}) - ({V_SPE})")
    print(f"This simplifies to the rational function p(d)/q(d):")
    
    # Print the final equation with each number explicitly
    p_coeffs = sympy.Poly(p_final, d).all_coeffs()
    q_coeffs = sympy.Poly(q_final, d).all_coeffs()
    
    p_str = f"({p_coeffs[0]})*d**2 + ({p_coeffs[1]})*d + ({p_coeffs[2]})"
    q_str = f"({q_coeffs[0]})*d**2 + ({q_coeffs[1]})*d + ({q_coeffs[2]})"
    print(f"p(d) = {p_final} which is {p_str}")
    print(f"q(d) = {q_final} which is {q_str}\n")
    
    # Step 6: Sum the absolute values of the coefficients
    sum_abs_coeffs_p = sum(abs(c) for c in p_coeffs)
    sum_abs_coeffs_q = sum(abs(c) for c in q_coeffs)
    total_sum = sum_abs_coeffs_p + sum_abs_coeffs_q
    
    print("--- Step 6: Sum of Absolute Coefficients ---")
    print(f"Coefficients of p(d) are: {p_coeffs}. Sum of absolute values: {sum_abs_coeffs_p}")
    print(f"Coefficients of q(d) are: {q_coeffs}. Sum of absolute values: {sum_abs_coeffs_q}")
    print(f"Total sum of absolute values of coefficients is {sum_abs_coeffs_p} + {sum_abs_coeffs_q} = {total_sum}")

if __name__ == '__main__':
    solve_bargaining_problem()