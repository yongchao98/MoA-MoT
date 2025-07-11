import numpy as np

def simulate_navier_stokes_2d():
    """
    Simulates the 2D incompressible Navier-Stokes equations using a pseudo-spectral method
    to demonstrate the decay of energy and enstrophy, illustrating the absence of finite-time blow-up.
    """
    # Simulation parameters
    N = 128          # Resolution (N x N grid)
    viscosity = 1.0  # Viscosity (nu), same as the Delta u term
    dt = 0.01        # Time step
    end_time = 5.0   # Simulation end time
    
    # Create the spatial grid
    L = 2 * np.pi
    dx = L / N
    x = np.arange(0, L, dx)
    X, Y = np.meshgrid(x, x, indexing='ij')

    # Create the wavenumber grid for Fourier space
    k = np.fft.fftfreq(N, d=dx) * L  # Wavenumbers scaled for [0, 2*pi) domain
    Kx, Ky = np.meshgrid(k, k, indexing='ij')
    K_sq = Kx**2 + Ky**2
    # Avoid division by zero at k=0
    K_sq_inv = np.zeros_like(K_sq)
    K_sq_inv[K_sq > 0] = 1.0 / K_sq[K_sq > 0]

    # Initial condition: a smooth vortex field
    # omega(x, y, t=0) = 2*cos(x)*sin(y)
    omega = 2 * np.cos(X) * np.sin(Y)
    omega_hat = np.fft.fft2(omega)

    # Time-stepping scheme uses an integrating factor for the linear (diffusion) term
    # and forward Euler for the nonlinear (advection) term.
    exp_factor = np.exp(-viscosity * K_sq * dt)
    
    # Main simulation loop
    time = 0.0
    step = 0
    print_interval = 10

    print("--- 2D Navier-Stokes Simulation ---")
    print("This simulation demonstrates that for smooth initial data, the solution does not blow up.")
    print("Instead, total kinetic energy and enstrophy decay over time.\n")
    print(f"{'Time':<10}{'Kinetic Energy':<20}{'Enstrophy':<20}")
    print("-" * 50)

    while time < end_time:
        if step % print_interval == 0:
            # --- Calculate diagnostics: energy and enstrophy ---
            
            # Get velocity from vorticity in Fourier space (Biot-Savart law)
            # u_hat = (i*ky/|k|^2, -i*kx/|k|^2) * omega_hat
            u_hat_x = 1j * Ky * K_sq_inv * omega_hat
            u_hat_y = -1j * Kx * K_sq_inv * omega_hat
            
            # Calculate total energy and enstrophy using Parseval's theorem
            # Note: The sum in Fourier space is proportional to the integral in real space.
            # Normalization factor (1/N^2)^2 from fft2, and (N^2) from sum. Total = 1/N^2.
            # We care about relative values, so precise scaling is not critical.
            energy = 0.5 * np.sum(np.abs(u_hat_x)**2 + np.abs(u_hat_y)**2) / N**2
            enstrophy = 0.5 * np.sum(np.abs(omega_hat)**2) / N**2
            
            print(f"{time:<10.2f}{energy:<20.6e}{enstrophy:<20.6e}")

        # --- Time-stepping ---
        # 1. Calculate nonlinear advection term (u . grad(omega)) in Fourier space
        
        # Get derivatives of omega in Fourier space
        omega_x_hat = 1j * Kx * omega_hat
        omega_y_hat = 1j * Ky * omega_hat
        
        # Transform fields to real space for pointwise multiplication
        u_x = np.fft.ifft2(u_hat_x)
        u_y = np.fft.ifft2(u_hat_y)
        omega_x = np.fft.ifft2(omega_x_hat)
        omega_y = np.fft.ifft2(omega_y_hat)
        
        # Calculate nonlinear term and transform back to Fourier space
        nonlinear_term = u_x * omega_x + u_y * omega_y
        nonlinear_hat = np.fft.fft2(nonlinear_term)
        
        # 2. Update vorticity in Fourier space
        omega_hat = exp_factor * (omega_hat - dt * nonlinear_hat)
        
        # 3. De-alias (optional, but good practice for stability)
        # This sets high-frequency modes to zero to remove aliasing errors.
        cutoff = N // 3
        omega_hat[np.abs(Kx) > cutoff] = 0
        omega_hat[np.abs(Ky) > cutoff] = 0

        time += dt
        step += 1

simulate_navier_stokes_2d()