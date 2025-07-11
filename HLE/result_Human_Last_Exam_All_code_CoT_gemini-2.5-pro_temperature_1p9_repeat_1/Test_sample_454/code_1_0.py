import math

def solve_lattice_gas_occupancy():
    """
    Calculates the average site occupancy for a lattice gas system
    using the mean-field approximation and prints the results.
    """
    # Step 1: Define problem parameters from the prompt.
    # We are given relationships, so we can work with dimensionless ratios.
    mu_over_kBT = 0.1
    epsilon_over_kBT = -1 / (2 * math.pi)
    z_horizontal = 4
    z_vertical = 8
    
    # In a mean-field model of a bulk system, the total coordination number
    # is the sum of its components.
    z = z_horizontal + z_vertical

    # Step 2: Display the mean-field self-consistency equation.
    print("In the mean-field approximation, the average site occupancy <n> is found by solving the self-consistency equation:")
    print("<n> = 1 / (1 + exp(-(μ - ε*z*<n>) / (k_B*T)))")
    print("\nThis equation can be rewritten with dimensionless parameters:")
    print("<n> = 1 / (1 + exp(-( (μ/(k_B*T)) - (ε/(k_B*T))*z*<n> )))")
    
    # Step 3: Print the equation with the specific numerical values.
    print("\nSubstituting the given values into the equation:")
    print(f"μ/(k_B*T) = {mu_over_kBT}")
    print(f"ε/(k_B*T) = -1/(2*π) = {epsilon_over_kBT:.5f}")
    print(f"z = z_horizontal + z_vertical = {z_horizontal} + {z_vertical} = {z}")

    # Display the final equation to be solved.
    print("\nThe specific equation to solve for <n> is:")
    print(f"<n> = 1 / (1 + exp(-({mu_over_kBT} - ({epsilon_over_kBT:.5f}) * {z} * <n>)))")
    print(f"<n> = 1 / (1 + exp(-({mu_over_kBT} + {abs(epsilon_over_kBT*z):.5f} * <n>)))")

    # Step 4: Solve the self-consistency equation numerically.
    # We use a fixed-point iteration method x_new = f(x).
    n = 0.5  # Initial guess for <n>
    tolerance = 1e-7 # Convergence tolerance
    max_iterations = 100 # Safety break for the loop

    for i in range(max_iterations):
        # The exponent in the self-consistency equation
        exponent = -(mu_over_kBT - epsilon_over_kBT * z * n)
        
        # Calculate the new value of n
        n_new = 1 / (1 + math.exp(exponent))
        
        # Check for convergence
        if abs(n_new - n) < tolerance:
            n = n_new
            break
        
        n = n_new
    else:
        print("\nWarning: Numerical solution did not converge within maximum iterations.")

    # Step 5: Print the final result, rounded to three decimal places.
    final_n_rounded = round(n, 3)

    print(f"\nAfter numerically solving the equation, the final average occupancy per site is:")
    print(f"{final_n_rounded}")
    
    # Return the answer in the specified format
    print(f"\n<<<{final_n_rounded}>>>")

# Run the calculation
solve_lattice_gas_occupancy()