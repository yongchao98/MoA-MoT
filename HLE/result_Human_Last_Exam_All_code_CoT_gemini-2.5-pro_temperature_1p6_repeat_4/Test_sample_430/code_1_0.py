import math

def calculate_adsorption():
    """
    Calculates the average number of adsorbed layers per site for a multi-layer
    lattice system using the Mean-Field Approximation and a set of simplifying
    assumptions based on the problem statement.
    """
    # --- Step 1: Define System Parameters from the problem ---
    # We are given energies and chemical potential in units of k_B*T.
    # We can work with beta*energy, where beta = 1/(k_B*T).
    beta_eps1 = 0.1
    beta_mu = 0.15
    z_ell = 4  # Lateral coordination number
    z_inter = 4 # Vertical interaction coordination number, used to determine eps_L
    k_max = 10 # Maximum number of layers. This value is assumed as it was not specified.
    
    # --- Step 2: Interpret Parameters and Make Assumptions ---
    
    # The lateral interaction energy is eps_l = (0.02)^k * k_B*T.
    # For k=k_max=10, eps_l is negligible, so we assume beta*eps_l = 0.
    beta_eps_l = 0.0 
    
    # The energies of subsequent layers (eps_2, eps_3, ...) are not given.
    # We assume they are all equal to a value eps_L derived from eps_1 and z_inter.
    # Assumption: eps_L = eps_1 / z_inter
    beta_eps_L = beta_eps1 / z_inter
    
    print("--- Model Parameters and Assumptions ---")
    print(f"Max layers (k_max): {k_max} (Assumed)")
    print(f"β*ε₁ = {beta_eps1}")
    print(f"β*μ = {beta_mu}")
    print(f"z_inter = {z_inter}")
    print(f"β*ε_L = β*ε₁/z_inter = {beta_eps_L} (Assumed for layers > 1)")
    print(f"β*ε_ℓ is assumed to be 0, so lateral interactions are ignored.\n")

    # --- Step 3: Calculate the x_i terms ---
    # These terms are the statistical weights for adding a particle to each layer.
    
    # x1 is for the first layer (adsorption on the surface)
    x1 = math.exp(beta_mu + beta_eps1)
    
    # xL is for subsequent layers (adsorption on another particle)
    xL = math.exp(beta_mu + beta_eps_L)

    print("--- Calculated Statistical Weights ---")
    print(f"x₁ = exp(β(μ+ε₁)) = exp({beta_mu + beta_eps1:.3f}) = {x1:.4f}")
    print(f"x_L = exp(β(μ+ε_L)) = exp({beta_mu + beta_eps_L:.3f}) = {xL:.4f}\n")

    # --- Step 4: Calculate the Single-Site Grand Partition Function (Ξ) and Numerator ---
    # We loop from j=1 to k_max to sum the terms.
    
    # The partition function starts with 1 (for the empty state, j=0)
    partition_function = 1.0
    
    # This sum will be the numerator in the expression for <k>
    numerator_sum = 0.0
    
    # This variable holds the term for the j-th state (x1 * xL^(j-1))
    # It is initialized to x1 for the j=1 case.
    term_j = x1
    
    # Calculate term for j=1 and add to sums
    partition_function += term_j
    numerator_sum += 1 * term_j
    
    # Loop for layers j = 2 to k_max
    for j in range(2, k_max + 1):
        # The term for layer j is the term for (j-1) multiplied by xL
        term_j *= xL
        partition_function += term_j
        numerator_sum += j * term_j
        
    print("--- Final Equation Components ---")
    print(f"Numerator Sum (Σ j * P(j) * Ξ): {numerator_sum:.4f}")
    print(f"Denominator (Grand Partition Function Ξ): {partition_function:.4f}\n")
    
    # --- Step 5: Calculate the Average Number of Layers (<k>) ---
    average_layers = numerator_sum / partition_function

    print("--- Final Result ---")
    print(f"The average number of adsorbed layers per site is <k> = Numerator / Ξ")
    print(f"<k> = {numerator_sum:.4f} / {partition_function:.4f} = {average_layers:.4f}")

    return average_layers

# Run the calculation and store the final answer
final_answer = calculate_adsorption()
print(f"\n<<<${final_answer:.4f}>>>")