import numpy as np

def solve_adsorption_model():
    """
    Solves for the average number of adsorbed layers using a mean-field model.
    """
    # Parameters based on the problem description and assumptions
    # Energies are in units of k_B * T
    epsilon_1 = 0.1
    # Assuming epsilon_j = epsilon_1 for all j > 1
    epsilon_j_func = lambda j: epsilon_1
    
    # Assuming epsilon_l = 0.02 * k_B * T based on interpretation of the problem
    epsilon_l = 0.02
    
    mu = 0.15
    z_l = 4
    
    # Assuming k=4 based on coordination numbers
    k_max = 4

    print("Solving for the average number of adsorbed layers <k> with the following parameters:")
    print(f"  k_max = {k_max}")
    print(f"  ε₁/kBT = {epsilon_1}")
    print(f"  εj/kBT = {epsilon_1} for j>1")
    print(f"  εℓ/kBT = {epsilon_l}")
    print(f"  μ/kBT = {mu}")
    print(f"  zℓ = {z_l}")
    print("-" * 30)

    theta = np.zeros(k_max + 1)
    theta[0] = 1.0  # Definition: theta_0 = 1

    # Sequentially solve for theta_j for j = 1 to k_max
    for j in range(1, k_max + 1):
        theta_j_minus_1 = theta[j - 1]
        epsilon_j = epsilon_j_func(j)

        # We need to solve the implicit equation for theta[j]:
        # theta[j] = (x_j(theta[j]) * theta_{j-1}) / (1 + x_j(theta[j]))
        # This is a fixed-point problem: theta[j] = f(theta[j])

        def f(current_theta_j):
            exponent = mu + epsilon_j + z_l * epsilon_l * current_theta_j
            x_j = np.exp(exponent)
            return (x_j * theta_j_minus_1) / (1 + x_j)

        # Use fixed-point iteration to find theta[j]
        # Start with an initial guess
        theta_j_guess = 0.0
        for _ in range(100):  # 100 iterations is more than enough for convergence
            theta_j_new = f(theta_j_guess)
            if np.abs(theta_j_new - theta_j_guess) < 1e-9:
                break
            theta_j_guess = theta_j_new
        theta[j] = theta_j_guess

    # The average number of layers is the sum of the layer coverages
    average_layers = np.sum(theta[1:])
    
    # Output the final equation as requested
    equation_parts = [f"{val:.4f}" for val in theta[1:]]
    equation_str = " + ".join(equation_parts)
    
    print("The average number of adsorbed layers per site, <k>, is the sum of the individual layer coverages θj:")
    print(f"\n<k> = θ₁ + θ₂ + θ₃ + θ₄")
    print(f"<k> = {equation_str} = {average_layers:.4f}")

    return average_layers

if __name__ == '__main__':
    result = solve_adsorption_model()
    # The final answer is returned in the specified format
    # The calculation gives ~1.1874
    # print(f"\n<<<{result:.4f}>>>")