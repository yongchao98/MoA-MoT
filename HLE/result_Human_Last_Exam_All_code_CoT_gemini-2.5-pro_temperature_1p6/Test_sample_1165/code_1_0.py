import numpy as np

def estimate_fluctuation_magnitude(epsilon=0.05, num_trials=500, num_x_points=1000):
    """
    Numerically estimates the maximum fluctuation magnitude R for the given ODE.

    Args:
        epsilon (float): The small parameter in the ODE.
        num_trials (int): The number of Monte Carlo simulations to run.
        num_x_points (int): The number of points to discretize the domain.

    Returns:
        float: The estimated value of R.
    """
    L = 1.0 / epsilon
    N = int(L) - 1
    
    # Create a grid for x values
    x_grid = np.linspace(0, L, num_x_points)
    
    # Array to store the solution y(x) for each trial
    y_samples = np.zeros((num_trials, num_x_points))

    # --- Green's function components ---
    # The Green's function G(x, s) solves L[G] = delta(x-s) with G(0,s)=G(L,s)=0
    # G(x, s) = c1(s) * (1 - exp(eps*x)) for x < s
    # G(x, s) = d2(s) * (exp(eps*x) - exp(eps*L)) for x > s
    e_val = np.exp(1.0)
    
    def c1(s_vals):
        # Note: np.exp(epsilon*s_vals) can be large, but L*epsilon=1, so it is bounded.
        numerator = np.exp(epsilon * s_vals) - np.exp(epsilon * L)
        denominator = epsilon * np.exp(epsilon * s_vals) * (1.0 - np.exp(epsilon*L))
        # Handle potential division by zero if s_vals are at boundaries, although unlikely
        denominator[denominator == 0] = 1e-9
        return numerator / denominator
        
    def d2(s_vals):
        numerator = 1.0 - np.exp(epsilon * s_vals)
        denominator = epsilon * np.exp(epsilon * s_vals) * (1.0 - np.exp(epsilon * L))
        # Handle potential division by zero
        denominator[denominator == 0] = 1e-9
        return numerator / denominator

    # --- Deterministic part of the solution y_h(x) ---
    y_h = (np.exp(epsilon * L) - np.exp(epsilon * x_grid)) / (np.exp(epsilon * L) - 1.0)

    print(f"Starting simulation for epsilon = {epsilon}")
    print(f"Domain L = {L}, Number of sources N = {N}")
    
    for i in range(num_trials):
        if (i + 1) % 100 == 0:
            print(f"  ... Running trial {i+1}/{num_trials}")
            
        # 1. Generate N ordered random points z_i in [0, L]
        z_i = np.sort(np.random.uniform(0, L, N))
        
        # 2. Calculate the stochastic part of the solution y_p(x) = eps^2 * sum(G(x, z_i))
        # We can vectorize this calculation over x_grid
        y_p = np.zeros_like(x_grid)
        
        # Calculate Green's function contributions for each z_i
        # For a given x, sum over contributions from all z_i
        c1_vals = c1(z_i)
        d2_vals = d2(z_i)
        
        for j, x in enumerate(x_grid):
            # z_i values where x < z_i
            mask_greater = z_i > x
            # z_i values where x > z_i
            mask_lesser = ~mask_greater
            
            # G(x, s) for x < s
            term1 = np.sum(c1_vals[mask_greater]) * (1 - np.exp(epsilon * x))
            
            # G(x, s) for x > s
            term2 = np.sum(d2_vals[mask_lesser]) * (np.exp(epsilon * x) - np.exp(epsilon * L))

            y_p[j] = term1 + term2
        
        y_p *= epsilon**2
        
        # 3. Full solution y(x) = y_h(x) + y_p(x)
        y_samples[i, :] = y_h + y_p

    # 4. Calculate variance and R
    # In the problem, y(0)=1 is a constant, so Var[y(x) - y(0)] = Var[y(x)]
    # As y_h is deterministic, Var[y(x)] = Var[y_p(x)]
    var_y = np.var(y_samples, axis=0)
    
    max_var = np.max(var_y)
    
    R = np.sqrt(max_var)
    
    # 5. Output final result
    print("\n--- Calculation Summary ---")
    print(f"Maximum variance found: max(Var[y(x)]) = {max_var:.6f}")
    print(f"The equation for R is: R = (max |Var[y(x)]|)^(1/2)")
    print(f"Plugging in the maximum variance:")
    print(f"R = ({max_var:.6f})^(1/2)")
    print(f"Final estimated magnitude of fluctuations R = {R:.6f}")
    
    # Theoretical prediction for scaling
    predicted_R_scaling = epsilon**0.5
    print(f"\nFor comparison, the theoretical scaling R ~ epsilon^0.5 gives a value proportional to {predicted_R_scaling:.6f}")


if __name__ == '__main__':
    # Set a small epsilon value
    # Smaller epsilon requires more memory/time (L is larger)
    epsilon_val = 0.04
    estimate_fluctuation_magnitude(epsilon=epsilon_val, num_trials=500, num_x_points=1000)
