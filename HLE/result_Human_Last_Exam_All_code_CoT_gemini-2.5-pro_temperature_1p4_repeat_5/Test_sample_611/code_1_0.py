import numpy as np

def solve():
    """
    Solves the problem by simulating the process and estimating the
    asymptotic sausage density p_infty.
    """

    def simulate_Xt(x0, dt, n_steps):
        """
        Simulates the process X_t using the Euler-Maruyama method.
        The SDE is dX_t = dW_t + (X_t / (|X_t|^2 * log|X_t|)) dt.
        """
        x = np.zeros(n_steps, dtype=np.complex128)
        x[0] = x0
        
        # Pre-generate random numbers for speed
        noise_real = np.random.normal(0, np.sqrt(dt / 2), n_steps - 1)
        noise_imag = np.random.normal(0, np.sqrt(dt / 2), n_steps - 1)
        
        for i in range(n_steps - 1):
            abs_x = np.abs(x[i])
            if abs_x <= 1.0:
                # Should not happen with the drift, but as a safeguard
                print(f"Warning: Process hit the unit disk at step {i}, stopping path.")
                x[i+1:] = x[i]
                break
            
            noise = noise_real[i] + 1j * noise_imag[i]
            drift = x[i] / (abs_x**2 * np.log(abs_x))
            x[i+1] = x[i] + noise + drift * dt
            
        return x

    def estimate_Vn(path, n, n_samples):
        """
        Estimates the volume fraction V_n using Monte Carlo integration.
        """
        center = n + 0j
        radius = n / 3.0
        
        # Generate random points uniformly in the disk B_n
        r_vals = radius * np.sqrt(np.random.rand(n_samples))
        theta_vals = 2 * np.pi * np.random.rand(n_samples)
        pts = center + r_vals * np.exp(1j * theta_vals)
        
        # Prune the path to the relevant section to speed up computation
        min_dist_from_origin = n - radius - 1
        max_dist_from_origin = n + radius + 1
        path_abs = np.abs(path)
        relevant_path_indices = np.where((path_abs > min_dist_from_origin) & (path_abs < max_dist_from_origin))[0]
        
        if len(relevant_path_indices) == 0:
            print(f"Warning: The path did not reach the vicinity of B_{n}.")
            return 0.0
            
        relevant_path = path[relevant_path_indices]

        # Use broadcasting for faster distance calculation in chunks
        covered_count = 0
        batch_size = 1000 
        for i in range(0, n_samples, batch_size):
            batch_pts = pts[i:i+batch_size]
            # distance matrix between batch_pts and relevant_path
            dist_matrix = np.abs(batch_pts[:, np.newaxis] - relevant_path)
            min_dists = np.min(dist_matrix, axis=1)
            covered_count += np.sum(min_dists <= 1.0)
            
        return covered_count / n_samples

    # --- Parameters ---
    np.random.seed(42) # for reproducibility
    x0 = 1.01 + 0j
    dt = 1e-4
    n_steps = 2 * 10**5 # T = 20
    n = 100
    n_samples_vn = 5 * 10**4 

    # --- Execution ---
    # print("Simulating the path of X_t...")
    path = simulate_Xt(x0, dt, n_steps)
    
    # print(f"Estimating V_n for n={n}...")
    p_inf_estimate = estimate_Vn(path, n, n_samples_vn)
    
    c = 2/3

    print(f"The asymptotic sausage density p_infty is estimated to be: {p_inf_estimate:.4f}")
    print(f"The threshold value is 2/3, which is approximately: {c:.4f}")

    if np.isclose(p_inf_estimate, c, atol=0.05):
        # This case suggests p_infty = 2/3, so the limiting probability is 1/2
        final_prob = 0.5
        print("The estimated density is close to 2/3. The limiting probability is 1/2.")
    elif p_inf_estimate > c:
        # p_infty > 2/3, so P(V_n > 2/3) -> 1
        final_prob = 1
        print("The estimated density is greater than 2/3. The limiting probability is 1.")
    else:
        # p_infty < 2/3, so P(V_n > 2/3) -> 0
        final_prob = 0
        print("The estimated density is less than 2/3. The limiting probability is 0.")
        
    print(f"\nFinal Equation: P(lim V_n > 2/3) = {final_prob}")

solve()