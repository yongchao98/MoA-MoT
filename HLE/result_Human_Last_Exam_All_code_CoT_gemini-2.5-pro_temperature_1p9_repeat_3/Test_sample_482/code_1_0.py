import sympy

def solve_critical_correlation():
    """
    Derives and prints the formula for the critical correlation 'C' required to
    balance potentiation and depression in the described neural network.
    """
    # Define symbolic variables to represent the parameters of the system.
    # C: The covariance between corresponding neurons in populations v and s, Cov(v_i, s_i). This is the 'correlation' we aim to find.
    # mu (μ): The average firing rate of the input neurons, <v_i> = <s_i> = μ.
    # theta (θ): The heterosynaptic offset constant in the plasticity rule.
    # sigma_sq (σ²): The variance of the input firing rates, Var(v_i) = Var(s_i) = σ².
    # W_v, W_s: The synaptic weights from the selective inputs v_i and s_i to the output neuron r_i.
    C, mu, theta, sigma_sq = sympy.symbols('C mu theta sigma_sq')
    W_v, W_s = sympy.symbols('W_v W_s', positive=True)

    # The activity of a selective output neuron 'i' is modeled as r_i = W_v*v_i + W_s*s_i.
    # The fixed-point conditions for the learning rules are:
    # 1. <r_i * (v_i - θ)> = 0
    # 2. <r_i * (s_i - θ)> = 0

    # We use the following statistical moments:
    # <v_i> = <s_i> = μ
    # <v_i²> = Var(v_i) + <v_i>² = σ² + μ²
    # <s_i²> = Var(s_i) + <s_i>² = σ² + μ²
    # <v_i * s_i> = Cov(v_i, s_i) + <v_i><s_i> = C + μ²

    # Expanding the first fixed-point equation:
    # <(W_v*v_i + W_s*s_i) * (v_i - θ)> = 0
    # W_v*<v_i²> + W_s*<s_i*v_i> - θ*W_v*<v_i> - θ*W_s*<s_i> = 0
    eq1 = W_v * (sigma_sq + mu**2) + W_s * (C + mu**2) - theta * W_v * mu - theta * W_s * mu

    # Expanding the second fixed-point equation:
    # <(W_v*v_i + W_s*s_i) * (s_i - θ)> = 0
    # W_v*<v_i*s_i> + W_s*<s_i²> - θ*W_v*<v_i> - θ*W_s*<s_i> = 0
    eq2 = W_v * (C + mu**2) + W_s * (sigma_sq + mu**2) - theta * W_v * mu - theta * W_s * mu

    # This is a homogeneous system of linear equations for W_v and W_s.
    # For a non-trivial solution to exist, the determinant of the coefficient matrix must be zero.
    coeff_Wv_eq1 = eq1.coeff(W_v)
    coeff_Ws_eq1 = eq1.coeff(W_s)
    coeff_Wv_eq2 = eq2.coeff(W_v)
    coeff_Ws_eq2 = eq2.coeff(W_s)

    determinant_eq = coeff_Wv_eq1 * coeff_Ws_eq2 - coeff_Ws_eq1 * coeff_Wv_eq2
    
    # Solving the determinant equation for C gives two possible conditions.
    # One corresponds to W_v = W_s (symmetric solution), and the other to W_v = -W_s (anti-symmetric solution).
    # We are interested in the symmetric solution where correlation drives selectivity.
    solutions_for_C = sympy.solve(determinant_eq, C)
    
    # Find the solution that does not equal sigma_sq, which corresponds to the symmetric W_v = W_s case.
    if solutions_for_C[0] == sigma_sq:
        symmetric_C_solution = solutions_for_C[1]
    else:
        symmetric_C_solution = solutions_for_C[0]
        
    # Format the result for clarity.
    final_solution = 2 * mu * (theta - mu) - sigma_sq

    # Create the final equation object C = ... for printing
    final_equation = sympy.Eq(C, final_solution)

    print("Based on the stability analysis of the learning rule at a fixed point, a specific relationship must hold between the input statistics.")
    print("This relationship defines the 'critical amount of correlation' (C) required to balance synaptic potentiation and depression.")
    print("\nThe critical amount of correlation is given by the following equation:")
    print("C = 2*μ*(θ - μ) - σ²")
    print("\nSymbolic representation:")
    sympy.pprint(final_equation)

if __name__ == '__main__':
    solve_critical_correlation()
    # For the final response, we are directly providing the code block.
    # The output from the code will be the formula derived.
    # So we don't need to capture and print <<<...>>> in this case, as the result is an equation.
    # Final answer is an equation, as calculated.
    C, mu, theta, sigma_sq = sympy.symbols('C mu theta sigma_sq')
    answer_expression = 2*mu*(theta-mu) - sigma_sq
    # The final deliverable needs to be a value. Since this problem is symbolic,
    # the symbolic expression itself is the value. I'll wrap it.
    final_answer = str(answer_expression)
    print(f"<<<{final_answer}>>>")