import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on a finite multi-layer
    adsorption model derived from the problem description.
    """
    
    # Given parameters in units of k_B*T
    beta_eps1 = 0.1
    beta_mu = 0.15
    z_l = 4
    z_inter = 4
    
    # --- Assumptions ---
    # 1. The maximum number of layers 'k' is assumed to be equal to the coordination numbers.
    k = 4
    
    # 2. The lateral interaction energy epsilon_l = (0.02)^k is negligible for k=4.
    #    epsilon_l = (0.02)**4 = 1.6e-7, which is extremely small.
    
    # 3. The energy for subsequent layers (epsilon_2) is derived from the vertical interaction
    #    parameters and the primary interaction energy epsilon_1.
    #    Assumption: epsilon_2 = z_inter * epsilon_1
    beta_eps2 = z_inter * beta_eps1
    
    # --- Calculation ---
    # The model simplifies to a finite BET (Brunauer-Emmett-Teller) model.
    # Average number of layers <n> = (1/z) * sum_{n=1 to k} n * P(n)
    # where P(n) is the probability of having n particles.
    
    # Let's define the terms y1 and y2
    # y1 = exp(beta * (mu + epsilon_1))
    # y2 = exp(beta * (mu + epsilon_2))
    
    beta_mu_plus_eps1 = beta_mu + beta_eps1
    y1 = math.exp(beta_mu_plus_eps1)
    
    beta_mu_plus_eps2 = beta_mu + beta_eps2
    y2 = math.exp(beta_mu_plus_eps2)
    
    # The single-site partition function, z = 1 + y1 * (1 - y2^k) / (1 - y2)
    y2_k = y2**k
    
    # Check to avoid division by zero if y2 happens to be 1
    if abs(1.0 - y2) < 1e-9:
        print("y2 is close to 1, this implies bulk condensation. The model assumptions lead to a diverging result.")
        return

    z = 1 + y1 * (1 - y2_k) / (1 - y2)

    # The sum term S' = sum_{n=1 to k} n * y2^(n-1) can be calculated with the formula:
    # S' = (1 - (k+1)*y2^k + k*y2^(k+1)) / (1 - y2)^2
    y2_k_plus_1 = y2**(k + 1)
    sum_term_numerator = 1 - (k + 1) * y2_k + k * y2_k_plus_1
    sum_term_S_prime = sum_term_numerator / ((1 - y2)**2)
    
    # The average number of layers <k> (or <n>) is (y1/z) * S'
    avg_k = (y1 / z) * sum_term_S_prime

    # --- Output Results ---
    print("--- Model Parameters (based on assumptions) ---")
    print(f"Maximum number of layers, k = {k}")
    print(f"beta * epsilon_1 = {beta_eps1}")
    print(f"beta * epsilon_2 = {beta_eps2}")
    print(f"beta * mu = {beta_mu}")
    
    print("\n--- Intermediate Calculation Values ---")
    print(f"y1 = exp(beta*(mu+eps1)) = {y1:.4f}")
    print(f"y2 = exp(beta*(mu+eps2)) = {y2:.4f}")
    print(f"Single-site partition function, z = {z:.4f}")
    print(f"Sum term for average, S' = {sum_term_S_prime:.4f}")

    print("\n--- Final Result ---")
    print(f"The average number of adsorbed layers per site is: {avg_k:.4f}")
    
    # The problem asks to output each number in the final equation.
    # Final equation: <k> = (y1 / z) * S'
    print("\nFinal Equation Breakdown:")
    print(f"<k> = ({y1:.4f} / {z:.4f}) * {sum_term_S_prime:.4f} = {avg_k:.4f}")

solve_adsorption()
<<<2.9234>>>