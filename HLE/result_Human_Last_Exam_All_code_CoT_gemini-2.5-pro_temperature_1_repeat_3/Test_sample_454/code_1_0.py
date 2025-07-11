import math

def solve_lattice_gas_occupancy():
    """
    Calculates the average occupancy per site for a lattice gas system
    using the mean-field approximation.
    """
    # --- Step 1: Define parameters from the problem statement ---
    # The temperature T = 300K is not needed explicitly as it cancels out.
    # We work with quantities relative to k_B*T.
    
    # Chemical potential over k_B*T
    mu_over_kBT = 0.1
    
    # Interaction energy over k_B*T
    epsilon_over_kBT = -1 / (2 * math.pi)
    
    # Coordination numbers
    z_horizontal = 4
    z_vertical = 8
    
    # The total coordination number z is the sum of horizontal and vertical neighbors
    z = z_horizontal + z_vertical
    
    # --- Step 2: Set up and solve the self-consistency equation numerically ---
    # The equation is: <n> = 1 / (exp(z * (epsilon/kBT) * <n> - mu/kBT) + 1)
    
    # Initial guess for the average occupancy <n>
    n = 0.5  
    
    # Iteration parameters
    tolerance = 1e-7
    max_iterations = 100
    
    for i in range(max_iterations):
        n_old = n
        
        # Calculate the exponent in the self-consistency equation
        exponent = z * epsilon_over_kBT * n_old - mu_over_kBT
        
        # Update the value of n
        n = 1 / (math.exp(exponent) + 1)
        
        # Check for convergence
        if abs(n - n_old) < tolerance:
            break
            
    # --- Step 3: Print the results ---
    
    print("Derivation Steps:")
    print("1. The self-consistency equation for average occupancy <n> in the mean-field approximation is:")
    print("   <n> = 1 / (exp(β * (z*ε*<n> - μ)) + 1)")
    print("\n2. Substituting β=1/(k_B*T), ε = -k_B*T/(2π), and μ = 0.1*k_B*T, the exponent becomes:")
    print(f"   β*(z*ε*<n> - μ) = z*({epsilon_over_kBT:.4f})*<n> - {mu_over_kBT}")
    
    print("\nFinal Equation:")
    print("The final equation to solve for <n> is:")
    # Print the equation with all numerical values included
    print(f"<n> = 1 / (exp({z} * ({epsilon_over_kBT:.4f}) * <n> - {mu_over_kBT}) + 1)")
    
    # Print the simplified form
    z_eps_term = z * epsilon_over_kBT
    print(f"<n> = 1 / (exp({z_eps_term:.4f} * <n> - {mu_over_kBT}) + 1)")
    
    print("\nResult:")
    # Print the final calculated value rounded to three decimal places
    print(f"The calculated average occupancy per site is: {n:.3f}")

# Execute the function
solve_lattice_gas_occupancy()
<<<0.848>>>