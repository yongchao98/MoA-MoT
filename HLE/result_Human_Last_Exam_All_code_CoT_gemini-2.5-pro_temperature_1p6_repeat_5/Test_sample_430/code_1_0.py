import math

def solve_monolayer_coverage():
    """
    Solves the self-consistent equation for monolayer coverage theta_1
    based on the Fowler-Guggenheim adsorption model.
    """
    # Parameters from the problem (in units of k_B*T)
    # k is set to 1 based on the stability analysis.
    k_max = 1
    epsilon_1 = 0.1
    mu = 0.15
    z_l = 4
    
    # Calculate lateral interaction energy for k=1
    epsilon_l = (0.02)**k_max

    # The equation to solve is:
    # theta = 1 / (1 + exp(z_l * epsilon_l * theta - (mu - epsilon_1)))
    # which is theta = 1 / (1 + exp(0.08 * theta - 0.05))
    
    # Print the equation we are solving
    c1 = z_l * epsilon_l
    c2 = mu - epsilon_1
    print("Solving the self-consistent Mean-Field equation for monolayer coverage theta:")
    print(f"theta = 1 / (1 + exp({c1:.3f} * theta - {c2:.3f}))")
    print(f"where:")
    print(f"  epsilon_1 / (k_B*T) = {epsilon_1}")
    print(f"  mu / (k_B*T) = {mu}")
    print(f"  z_l = {z_l}")
    print(f"  epsilon_l / (k_B*T) = (0.02)^k = {epsilon_l} (for k={k_max})")
    print("-" * 30)

    # Iterative solver for theta_1
    theta = 0.5  # Initial guess
    for i in range(100): # Iterate to find the fixed point
        theta_new = 1 / (1 + math.exp(z_l * epsilon_l * theta - (mu - epsilon_1)))
        # Check for convergence
        if abs(theta_new - theta) < 1e-9:
            break
        theta = theta_new

    # The average number of adsorbed layers is the coverage theta_1
    avg_k = theta
    
    print(f"The calculated average number of adsorbed layers per site is:")
    print(f"<k> = theta_1 = {avg_k}")
    
    return avg_k

if __name__ == "__main__":
    final_answer = solve_monolayer_coverage()
    # The final numerical answer in the required format
    print(f"\n<<<__{final_answer:.4f}__>>>")
