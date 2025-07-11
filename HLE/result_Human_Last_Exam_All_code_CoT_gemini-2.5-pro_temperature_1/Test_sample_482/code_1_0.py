import sympy

def solve_critical_correlation():
    """
    This function derives the critical correlation required to balance potentiation
    and depression in the described neural network.
    """
    # 1. Define the symbolic variables from the problem description.
    # N: Number of neurons in each input layer (assuming N_v = N_s = N)
    # mu: Average activation rate of input neurons (<v_k> = <s_k> = mu)
    # theta: Heterosynaptic offset constant in the plasticity rule
    # c_vs: Covariance between corresponding input neurons v_k and s_k
    N, mu, theta, c_vs = sympy.symbols('N mu theta c_vs')
    rho_vs = sympy.Symbol('rho_vs')

    # 2. Explain the stability condition.
    print("The balance between potentiation and depression is achieved when the total synaptic")
    print("weight, W_tot, onto an output neuron is stable (i.e., its growth rate is zero).\n")

    # 3. State the growth rate of the total synaptic weight.
    # This is derived from the system's dynamic equations for the weights. The derivation
    # shows that d(W_tot)/dt is proportional to W_tot, with a growth rate lambda_sum.
    # lambda_sum = <r_i * (v_k - theta)> + <r_i * (s_k - theta)> averaged over all synapses
    # After calculation, this growth rate is found to be:
    lambda_sum = N * mu**2 - N * theta * mu + mu - mu**2 + c_vs
    
    print("The growth rate (lambda_sum) of the total synaptic weight is given by:")
    print(f"lambda_sum = {lambda_sum}\n")

    # 4. Solve for the critical covariance by setting the growth rate to zero.
    print("To achieve stability, we set this growth rate to zero and solve for the critical covariance (c_vs_crit).\n")
    
    stability_equation = sympy.Eq(lambda_sum, 0)
    critical_c_vs = sympy.solve(stability_equation, c_vs)[0]
    
    print("The critical covariance is:")
    print(f"c_vs_crit = {critical_c_vs}\n")
    
    # 5. Convert the critical covariance to a critical correlation.
    # The correlation is rho = Cov(v,s) / (sigma_v * sigma_s).
    # For the described input process, the variance is mu * (1 - mu).
    input_variance = mu * (1 - mu)
    
    print("The correlation coefficient (rho_vs) is the covariance divided by the input variance.")
    print(f"Input Variance = {input_variance}\n")
    
    critical_rho_vs = sympy.simplify(critical_c_vs / input_variance)

    # 6. Display the final equation for the critical correlation.
    print("The final equation for the critical amount of correlation is:")
    print("----------------------------------------------------------")
    # We construct the final equation string manually for clear output
    final_eq_str = f"{rho_vs} = (N * ({theta} - {mu})) / (1 - {mu}) - 1"
    print(final_eq_str)
    print("----------------------------------------------------------\n")
    
    print("Where:")
    print(f"  {rho_vs}:\tThe critical correlation coefficient between inputs v and s.")
    print(f"  {N}:\tThe number of neurons in an input layer.")
    print(f"  {mu}:\tThe average activation rate of input neurons.")
    print(f"  {theta}:\tThe heterosynaptic offset constant.")

solve_critical_correlation()