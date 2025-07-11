import sympy

def display_critical_correlation_analysis():
    """
    This function explains and calculates the critical correlation needed to balance
    potentiation and depression in the described neural network, allowing for
    the development of input selectivity.
    """

    # --- 1. Define symbolic parameters for clarity ---
    # Using sympy for symbolic representation to make the equations clear.
    N = sympy.Symbol('N')  # Number of neurons per input layer (N_v = N_s = N)
    sigma_sq = sympy.Symbol('σ^2')  # Variance of input neuron activity
    C_crit = sympy.Symbol('C_crit')  # Critical covariance
    rho_crit = sympy.Symbol('ρ_crit')  # Critical correlation coefficient

    # --- 2. Provide theoretical background ---
    print("Theoretical Derivation of the Critical Correlation")
    print("-" * 50)
    print("The stability of weights W^v and W^s requires the expected change to be zero.")
    print("From the learning rules, this leads to a fixed-point condition:")
    print("  E[r_i * (v_k - θ)] = 0   and   E[r_i * (s_k - θ)] = 0")
    print("This implies that E[r_i * v_k] = E[r_i * s_k], which simplifies to E[r_i * (v_k - s_k)] = 0.")
    print("\nBy substituting r_i and applying mean-field analysis, we arrive at the bifurcation equation:")
    print("  (<W^v> - <W^s>) * (σ^2 - N * C) = 0")
    print("\nFor selectivity to develop, the average weights must differ (<W^v> != <W^s>).")
    print("This is only possible if the second term is zero, which defines the critical point.")
    print(f"  {sigma_sq} - {N} * C = 0")
    print("\nWhere:")
    print(f"  C: Covariance between input populations v and s, cov(v, s)")
    print(f"  {sigma_sq}: Variance of the input activities, var(v) = var(s)")
    print(f"  {N}: Number of neurons in each input population (assuming N_v = N_s = N)")
    print("-" * 50)

    # --- 3. Define and print the final equations with example values ---
    
    # Use example numerical values for calculation
    N_val = 100
    sigma_sq_val = 10.0  # Example: Variance for a process with a mean rate of 10 Hz

    print("\nCalculating the Critical Covariance (C_crit)")
    # The final equation for critical covariance
    critical_covariance_eq = sympy.Eq(C_crit, sigma_sq / N)
    print("The critical covariance is given by the equation:")
    print(f"  {sympy.pretty(critical_covariance_eq, use_unicode=False)}")

    # Calculate the numerical result
    c_crit_val = sigma_sq_val / N_val
    print("\nUsing example values:")
    print(f"  {C_crit} = {sigma_sq_val} / {N_val}")
    print(f"Result: {C_crit} = {c_crit_val}")
    print("-" * 50)

    print("\nCalculating the Critical Correlation Coefficient (ρ_crit)")
    # The final equation for the critical correlation coefficient
    # rho = C / sigma_sq
    critical_rho_eq = sympy.Eq(rho_crit, (sigma_sq / N) / sigma_sq)
    simplified_rho_eq = sympy.Eq(rho_crit, 1 / N)

    print("The critical correlation coefficient (ρ = C / σ²) is given by the equation:")
    print(f"  {sympy.pretty(critical_rho_eq, use_unicode=False)}")
    print("Which simplifies to:")
    print(f"  {sympy.pretty(simplified_rho_eq, use_unicode=False)}")
    
    # Calculate the numerical result
    rho_crit_val = 1 / N_val
    print("\nUsing example values:")
    print(f"  {rho_crit} = 1 / {N_val}")
    print(f"Result: {rho_crit} = {rho_crit_val}")
    print("-" * 50)

# Execute the analysis
display_critical_correlation_analysis()