import math

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site based on a
    multilayer adsorption model with specific simplifying assumptions.
    """
    # Dimensionless parameters from the problem statement (in units of k_B*T)
    # Based on the plan, we assume mu is negative to ensure physical stability.
    beta_mu = -0.15
    beta_eps1 = 0.1
    
    # As per the plan, for k -> infinity, beta_eps_ell and beta_eps_i (for i>1) are zero.
    
    print("Step 1: Calculate the equilibrium constant K1 for the first layer.")
    # K1 = exp(β(μ + ε₁))
    K1 = math.exp(beta_mu + beta_eps1)
    print(f"βμ = {beta_mu}, βε₁ = {beta_eps1}")
    print(f"K1 = exp({beta_mu} + {beta_eps1}) = {K1:.4f}\n")
    
    print("Step 2: Calculate the coverage of the first layer, θ₁.")
    # θ₁ / (1 - θ₁) = K1  =>  θ₁ = K1 / (1 + K1)
    theta1 = K1 / (1 + K1)
    print(f"θ₁ = K1 / (1 + K1) = {K1:.4f} / (1 + {K1:.4f}) = {theta1:.4f}\n")

    print("Step 3: Calculate the equilibrium constant K_upper for subsequent layers.")
    # K_upper = exp(β(μ + ε_i)) where ε_i = 0 for i >= 2
    K_upper = math.exp(beta_mu)
    print(f"K_upper = exp({beta_mu}) = {K_upper:.4f}\n")

    print("Step 4: Calculate the ratio 'r' for the geometric series of layer coverages.")
    # r = θ_i / θ_{i-1} = K_upper / (1 + K_upper)
    r = K_upper / (1 + K_upper)
    print(f"r = K_upper / (1 + K_upper) = {K_upper:.4f} / (1 + {K_upper:.4f}) = {r:.4f}\n")

    print("Step 5: Calculate the average number of layers <n> by summing the geometric series.")
    # <n> = θ₁ * (1 + r + r^2 + ...) = θ₁ / (1 - r)
    average_n = theta1 / (1 - r)
    print(f"<n> = θ₁ / (1 - r) = {theta1:.4f} / (1 - {r:.4f}) = {theta1:.4f} / {1-r:.4f}\n")
    
    print(f"The final calculated average number of adsorbed layers per site is: {average_n:.4f}")
    
    # Return the final value for the answer tag
    return average_n

final_answer = calculate_average_layers()
# print(f'<<<{final_answer:.4f}>>>')