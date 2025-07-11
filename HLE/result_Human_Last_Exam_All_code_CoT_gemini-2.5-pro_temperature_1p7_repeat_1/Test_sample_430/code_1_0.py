import math

def solve_adsorption_model():
    """
    Solves the self-consistent mean-field model for multilayer adsorption.
    """
    # Parameters from the problem description, converted to dimensionless form
    # (i.e., in units of k_B*T).
    epsilon_1_prime = 0.1
    mu_prime = 0.15
    z_l = 4

    # Assumption: The term (0.02)^k for epsilon_l is a typo for a constant value.
    epsilon_l_prime = 0.02

    # Assumption: Maximum number of layers for calculation, chosen to be large
    # enough for the sums to converge.
    K_max = 20

    # Self-consistency loop to find <k> (denoted as theta)
    # Let theta = <k>
    # The equation to solve is theta = f(theta)
    
    # Pre-calculate constant terms
    c1 = epsilon_1_prime - mu_prime
    c2 = z_l * epsilon_l_prime

    # Iterative solver
    theta = 1.0  # Initial guess for <k>
    for _ in range(200):  # Iterate up to 200 times
        
        # Calculate terms for the sums based on the current theta
        exp_terms = [math.exp(-(c1 + k * c2 * theta)) for k in range(1, K_max + 1)]

        # Calculate the numerator for the <k> expression
        numerator = sum(k * exp_terms[k - 1] for k in range(1, K_max + 1))

        # Calculate the single-site partition function z_s
        z_s = 1.0 + sum(exp_terms)

        # Avoid division by zero
        if z_s == 0:
            theta_new = 0.0
        else:
            # Calculate the new value of theta
            theta_new = numerator / z_s
        
        # Check for convergence
        if abs(theta_new - theta) < 1e-9:
            theta = theta_new
            break
        
        theta = theta_new

    # Final converged value of <k>
    final_avg_k = theta

    print(f"Based on the mean-field model derived under the specified assumptions:")
    print(f"The average number of adsorbed layers per site is calculated to be:\n")
    print(f"<k> = {final_avg_k:.6f}\n")
    
    # As requested, output the final numbers in the equation for the last iteration
    final_exp_terms = [math.exp(-(c1 + k * c2 * final_avg_k)) for k in range(1, K_max + 1)]
    final_z_s = 1.0 + sum(final_exp_terms)
    final_numerator = sum(k * final_exp_terms[k - 1] for k in range(1, K_max + 1))
    
    print("The final self-consistent equation with the calculated numbers is:")
    print(f"<{K_max}> = ({final_numerator:.6f}) / ({final_z_s:.6f})")
    
    print("\nFor clarity, here are the first few terms of the numerator sum:")
    print(f"Numerator = 1 * exp(-({c1:.2f} + 1*{c2:.2f}*{final_avg_k:.4f})) + 2 * exp(...) + ...")
    term1_val = final_exp_terms[0]
    term2_val = final_exp_terms[1]
    print(f"          = 1 * {term1_val:.4f} + 2 * {term2_val:.4f} + ...")

solve_adsorption_model()