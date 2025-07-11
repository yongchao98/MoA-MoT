import math

def solve_lattice_gas():
    """
    Calculates the average occupancy per site for a lattice gas model
    using the mean-field approximation.
    """
    # Given parameters
    mu_over_kBT = 0.1
    z_horizontial = 4
    z_vertical = 8
    
    # Step 1: Calculate the total coordination number z
    z = z_horizontial + z_vertical
    
    # Step 2: Define the dimensionless parameters for the self-consistent equation.
    # The equation is: <n> = 1 / (1 + exp(-beta * (mu - z * epsilon * <n>)))
    # Let A = beta * mu and B = beta * z * epsilon
    
    # A = beta * mu = (1/(kB*T)) * (0.1 * kB*T)
    A = mu_over_kBT
    
    # B = beta * z * epsilon = (1/(kB*T)) * z * (-kB*T / (2*pi))
    B = z * (-1 / (2 * math.pi))
    
    # Step 3: Solve the self-consistent equation using fixed-point iteration.
    # The equation is: n_avg = 1 / (1 + exp(-(A - B * n_avg)))
    n_avg = 0.5  # Initial guess for the average occupancy
    tolerance = 1e-7
    
    for i in range(100):  # Iterate up to 100 times to find the fixed point
        n_new = 1 / (1 + math.exp(-(A - B * n_avg)))
        if abs(n_new - n_avg) < tolerance:
            n_avg = n_new
            break
        n_avg = n_new
        
    # Step 4: Output the explanation, the final equation, and the result.
    print("In the mean-field approximation, the average occupancy <n> is found by solving the self-consistent equation:")
    print("<n> = 1 / (1 + exp(-beta * (mu - z * epsilon * <n>)))")
    print("\nUsing the given parameters, the final equation to solve is:")
    print(f"<n> = 1 / (1 + exp(-({A:.3f} - ({B:.3f}) * <n>)))")
    
    print(f"\nWhere:")
    print(f"  beta * mu = {A:.3f}")
    print(f"  z = {z}")
    print(f"  beta * z * epsilon = {B:.3f}")

    print(f"\nSolving this equation numerically yields the average occupancy per site:")
    print(f"<n> = {n_avg:.3f}")

solve_lattice_gas()
<<<0.848>>>