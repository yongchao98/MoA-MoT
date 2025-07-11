import math

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site using the
    finite-layer BET model under a set of assumptions for the missing parameters.
    """
    # Given parameters in units of k_B*T
    beta_eps1 = 0.1
    beta_mu = 0.15

    # Assumptions for the missing parameters k and epsilon_s
    # Assume k is equal to the given vertical coordination number z_inter
    k = 4
    # Assume the subsequent layer interaction energy epsilon_s is equal to the chemical potential mu
    beta_eps_s = beta_mu

    # The lateral interaction energy epsilon_l is negligible for k>=2
    beta_eps_l = (0.02)**k
    z_l = 4
    
    # Calculate the parameters for the finite-BET equation
    # c is related to the difference in adsorption energy between the first and subsequent layers
    c = math.exp(beta_eps1 - beta_eps_s)
    # xs is related to the chemical potential and the adsorption energy of subsequent layers
    xs = math.exp(beta_eps_s + beta_mu)

    # Calculate the average number of layers <h> using the finite-layer BET formula:
    # <h> = (c*xs/(1-xs)) * (1 - (k+1)*xs^k + k*xs^(k+1)) / (1 + (c-1)*xs - c*xs^(k+1))
    
    xs_k = xs**k
    xs_k_plus_1 = xs**(k+1)

    prefactor = (c * xs) / (1 - xs)
    
    numerator = 1 - (k + 1) * xs_k + k * xs_k_plus_1
    
    denominator = 1 + (c - 1) * xs - c * xs_k_plus_1
    
    # Handle the case where the denominator might be zero, though unlikely with these parameters
    if denominator == 0:
        print("Error: Division by zero in the BET formula.")
        return

    h_avg = prefactor * (numerator / denominator)

    # Print the parameters used in the final equation
    print(f"--- Parameters for Calculation ---")
    print(f"Adsorption energy for the first layer (in units of k_B*T): {beta_eps1}")
    print(f"Adsorption energy for subsequent layers (in units of k_B*T): {beta_eps_s:.4f} (Assumed)")
    print(f"Chemical potential (in units of k_B*T): {beta_mu}")
    print(f"Maximum number of layers, k: {k} (Assumed)")
    print(f"Lateral interaction term z_l*epsilon_l (in k_B*T): {z_l*beta_eps_l:.2e} (Negligible)")
    print("\n--- Intermediate Values for BET Equation ---")
    print(f"c = exp(beta*(epsilon_1 - epsilon_s)) = {c:.4f}")
    print(f"x_s = exp(beta*(epsilon_s + mu)) = {xs:.4f}")
    
    # Print the final result in the requested format
    print("\n--- Final Result ---")
    print(f"The average number of adsorbed layers per site is: {h_avg:.4f}")
    print(f"\n<<< {h_avg:.4f} >>>")


calculate_average_layers()