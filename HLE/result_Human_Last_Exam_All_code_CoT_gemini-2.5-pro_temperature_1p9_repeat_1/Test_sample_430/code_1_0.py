import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site
    based on a simplified multi-layer adsorption model.
    """
    # Given parameters normalized by k_B*T
    # mu / (k_B*T)
    mu_norm = 0.15
    # epsilon_1 / (k_B*T)
    eps1_norm = 0.1
    
    # Assumptions based on the problem statement
    # k_max is assumed to be 4, from z_inter = 4
    k_max = 4
    # epsilon_l is negligible due to the (0.02)^k factor
    # epsilon_i for i>=2 is assumed to be 0 as epsilon_inter is not given.
    
    # Calculate the statistical weights x_i for each layer
    # x_i = exp(beta * (mu + epsilon_i))
    x1 = math.exp(mu_norm + eps1_norm)
    x2 = math.exp(mu_norm)

    # Calculate the partition function z for a single site
    # z = 1 + x1 + x1*x2 + x1*x2^2 + ...
    weights = [1.0] # a_0 term
    current_weight = 1.0
    for h in range(1, k_max + 1):
        if h == 1:
            current_weight *= x1
        else:
            current_weight *= x2
        weights.append(current_weight)

    z = sum(weights)
    
    # Calculate the average number of layers <k>
    # <k> = (1/z) * sum_{h=1 to k_max} h * weight_h
    h_avg_numerator = 0
    for h in range(1, k_max + 1):
        h_avg_numerator += h * weights[h]
        
    h_avg = h_avg_numerator / z

    # Output the steps and the final equation with numerical values
    print("Based on the assumptions (k_max=4, ε_l=0, ε_i=0 for i>1), we calculate the average number of layers <k>.\n")
    print(f"The equation for the average number of layers is:")
    print(f"<k> = (Σ_{{h=1 to {k_max}}} h * w_h) / (Σ_{{h=0 to {k_max}}} w_h)\n")
    
    print("Component values:")
    print(f"k_max = {k_max}")
    print(f"βμ = {mu_norm}")
    print(f"βε₁ = {eps1_norm}")
    print(f"x₁ = exp(βμ + βε₁) = exp({mu_norm} + {eps1_norm:.1f}) = {x1:.4f}")
    print(f"x₂ = exp(βμ) = exp({mu_norm}) = {x2:.4f}\n")
    
    # Display weights
    print("Statistical weights for h layers (w_h):")
    for h in range(k_max + 1):
        print(f"w_{h} = {weights[h]:.4f}")

    # Numerator calculation
    num_str = " + ".join([f"{h}*{weights[h]:.4f}" for h in range(1, k_max + 1)])
    print(f"\nNumerator = Σ h*w_h = {num_str} = {h_avg_numerator:.4f}")
    
    # Denominator (z) calculation
    den_str = " + ".join([f"{w:.4f}" for w in weights])
    print(f"Denominator z = Σ w_h = {den_str} = {z:.4f}")

    # Final result
    print(f"\nAverage number of layers <k> = Numerator / Denominator = {h_avg_numerator:.4f} / {z:.4f}")
    print(f"Final Answer: {h_avg:.4f}")
    print(f"\n<<<${h_avg:.4f}>>>")

solve_adsorption()