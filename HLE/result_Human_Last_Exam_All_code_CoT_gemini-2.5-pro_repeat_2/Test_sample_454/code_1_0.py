import numpy as np

def calculate_average_occupancy():
    """
    This script calculates the average occupancy per site for a lattice gas
    using the mean-field approximation.
    """
    
    # Step 1: Explain the model and derive the self-consistency equation.
    print("Step 1: Deriving the Self-Consistency Equation")
    print("--------------------------------------------------")
    print("We model the system as a lattice gas where each site can be empty or occupied.")
    print("In the mean-field approximation, the interaction of a particle with its 'z' neighbors (energy epsilon) is replaced by an average interaction field.")
    print("The effective energy for a particle on a site, relative to the chemical potential mu, is:")
    print("E_eff = z * epsilon * <n> - mu")
    print("where <n> is the average occupancy per site we want to find.")
    
    print("\nThe probability of a site being occupied is given by the Fermi-Dirac-like distribution:")
    print("<n> = 1 / (exp(beta * E_eff) + 1), where beta = 1/(k_B*T).")
    
    print("\nSubstituting E_eff, we get the self-consistency equation for <n>:")
    print("<n> = 1 / (exp(beta * (z * epsilon * <n> - mu)) + 1)")
    print("")

    # Step 2: Define and substitute the given parameters.
    print("Step 2: Substituting Parameters into the Equation")
    print("--------------------------------------------------")
    
    z_horizontial = 4
    z_vertical = 8
    # The total coordination number z is the sum of horizontal and vertical neighbors.
    z = z_horizontial + z_vertical
    
    # The parameters are given in units of k_B*T. We can work with the dimensionless product beta*E.
    # beta * mu = (1/(k_B*T)) * (0.1 * k_B*T)
    beta_mu = 0.1
    # beta * epsilon = (1/(k_B*T)) * (-k_B*T / (2*pi))
    beta_epsilon = -1 / (2 * np.pi)
    
    # Calculate the coefficient for <n> in the exponent
    beta_z_epsilon = z * beta_epsilon

    print(f"Total coordination number, z = z_horizontial + z_vertical = {z_horizontial} + {z_vertical} = {z}")
    print(f"The dimensionless chemical potential, beta * mu = {beta_mu}")
    print(f"The dimensionless interaction energy, beta * epsilon = {beta_epsilon:.4f}")
    
    print("\nThe final numerical equation to solve is:")
    print(f"<n> = 1 / (exp({beta_z_epsilon:.5f} * <n> - {beta_mu}) + 1)")
    print("")

    # Step 3: Solve the equation numerically using an iterative method.
    print("Step 3: Solving the Equation Numerically")
    print("--------------------------------------------------")
    
    n_guess = 0.5  # Initial guess for <n>
    tolerance = 1e-7
    max_iterations = 100
    
    print(f"Starting iterative solution with initial guess <n> = {n_guess}")
    for i in range(max_iterations):
        exponent = beta_z_epsilon * n_guess - beta_mu
        n_new = 1 / (np.exp(exponent) + 1)
        
        if abs(n_new - n_guess) < tolerance:
            print(f"Converged after {i+1} iterations.")
            n_guess = n_new
            break
            
        n_guess = n_new
    
    final_n = n_guess
    
    print(f"\nThe calculated average occupancy per site is <n> = {final_n:.5f}")
    print("\nRounding to three decimal places, the final answer is:")
    print(f"{final_n:.3f}")
    
    return final_n

if __name__ == '__main__':
    result = calculate_average_occupancy()
    print(f"\n<<<0.848>>>")
