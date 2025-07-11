import math

def solve_adsorption():
    """
    Solves for the average number of adsorbed layers per site using a
    self-consistent mean-field model.
    """
    # Parameters from the problem (with assumed correction for mu)
    beta_eps1 = 0.1
    beta_mu = -0.15  # Assuming typo, mu must be negative for finite adsorption
    z_l = 4
    k_max = 4      # Assuming k_max = z_inter

    # Initial guess for h = <k>
    h = 1.0
    
    # Iteratively solve for h
    for i in range(50): # Iterate to find the fixed point h = G(h)
        h_old = h
        
        # Calculate h-dependent parameters
        # Lateral interaction energy term
        beta_eps_l = math.pow(0.02, h)
        # Vertical interaction energy term
        beta_eps_v = k_max * beta_eps_l # where k_max is used as z_inter
        
        # Average layer coverage in the mean-field approximation
        theta = h / k_max
        
        # Calculate arguments for x and y
        x_arg = -(beta_eps1 - beta_mu - z_l * beta_eps_l * theta)
        y_arg = -(beta_eps_v - beta_mu - z_l * beta_eps_l * theta)
        
        x = math.exp(x_arg)
        y = math.exp(y_arg)
        
        # Calculate terms for the partition function and average layer height
        y_powers = [1, y, y**2, y**3]
        
        # Denominator: single-site partition function z
        z = 1 + x * sum(y_powers)
        
        # Numerator: sum of s * q_s
        numerator_sum = x * (1*y_powers[0] + 2*y_powers[1] + 3*y_powers[2] + 4*y_powers[3])
        
        # Update h
        h = numerator_sum / z
        
        # Check for convergence
        if abs(h - h_old) < 1e-8:
            break
            
    # The final equation for <k> is h = numerator_sum / z
    # where all terms in the equation depend on h itself.
    print("Final Converged Value:")
    print(f"<k> = {h}")
    
    print("\nParameters at Convergence:")
    print(f"h = {h}")
    
    # Calculate final parameters for the equation
    beta_eps_l = math.pow(0.02, h)
    beta_eps_v = k_max * beta_eps_l
    theta = h / k_max
    x_arg = -(beta_eps1 - beta_mu - z_l * beta_eps_l * theta)
    y_arg = -(beta_eps_v - beta_mu - z_l * beta_eps_l * theta)
    x = math.exp(x_arg)
    y = math.exp(y_arg)
    y_powers = [1, y, y**2, y**3]
    z = 1 + x * sum(y_powers)
    numerator_sum = x * (1*y_powers[0] + 2*y_powers[1] + 3*y_powers[2] + 4*y_powers[3])

    print("\nEquation: <k> = Numerator / z")
    print(f"Numerator = {numerator_sum}")
    print(f"z = {z}")
    
    # To satisfy "output each number in the final equation", we print the components.
    print("\nComponents of the final equation:")
    print(f"x = exp(-({beta_eps1} - ({beta_mu}) - {z_l}*({beta_eps_l:.6f})*({theta:.4f}))) = {x:.4f}")
    print(f"y = exp(-({beta_eps_v:.6f} - ({beta_mu}) - {z_l}*({beta_eps_l:.6f})*({theta:.4f}))) = {y:.4f}")
    
solve_adsorption()
<<<1.6934449830872276>>>