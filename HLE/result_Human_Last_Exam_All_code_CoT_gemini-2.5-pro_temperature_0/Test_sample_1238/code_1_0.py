import numpy as np
from scipy.fft import fft2, ifft2, fftfreq

def solve_and_observe_2d_navier_stokes():
    """
    Simulates the 2D incompressible Navier-Stokes equation to demonstrate
    the absence of finite-time blow-up.

    This function solves the vorticity equation:
    ∂_t ω + u ⋅ ∇ω = ν Δω
    where u is the velocity and ω = ∇ × u is the vorticity.
    A pseudo-spectral method is used on a periodic domain [0, 2π] x [0, 2π].
    """
    # --- 1. Simulation Parameters ---
    N = 128          # Resolution (N x N grid)
    nu = 1.0         # Viscosity (as in the problem statement, Δu)
    dt = 0.01        # Time step
    end_time = 5.0   # Total simulation time
    
    # --- 2. Grid and Wavenumbers ---
    L = 2 * np.pi    # Domain size
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Wavenumbers
    k_vals = fftfreq(N, d=L/(2*np.pi))
    Kx, Ky = np.meshgrid(k_vals, k_vals, indexing='ij')
    K2 = Kx**2 + Ky**2
    # Avoid division by zero at k=0
    K2_inv = np.zeros_like(K2)
    K2_inv[K2 > 0] = 1.0 / K2[K2 > 0]

    # De-aliasing mask (2/3 rule)
    k_max = np.max(k_vals)
    dealias_mask = (np.abs(Kx) < (2./3.)*k_max) & (np.abs(Ky) < (2./3.)*k_max)

    # --- 3. Initial Condition (Taylor-Green vortex) ---
    # u_0 = (cos(x)sin(y), -sin(x)cos(y))
    # ω_0 = ∂_x u_2 - ∂_y u_1 = -cos(x)cos(y) - cos(x)cos(y) = -2cos(x)cos(y)
    omega_phys = -2.0 * np.cos(X) * np.cos(Y)
    omega_hat = fft2(omega_phys) # Vorticity in Fourier space

    # --- 4. RHS of the ODE in Fourier space ---
    def compute_rhs(w_hat):
        """Computes the RHS of the vorticity equation in Fourier space."""
        # Get velocity from vorticity in Fourier space (Biot-Savart law)
        # u_hat = (i*ky/|k|^2 * w_hat, -i*kx/|k|^2 * w_hat)
        u_hat = 1j * Ky * K2_inv * w_hat
        v_hat = -1j * Kx * K2_inv * w_hat

        # Derivatives of vorticity in Fourier space
        omega_x_hat = 1j * Kx * w_hat
        omega_y_hat = 1j * Ky * w_hat

        # Transform to physical space for nonlinear term calculation
        u_phys = ifft2(u_hat).real
        v_phys = ifft2(v_hat).real
        omega_x_phys = ifft2(omega_x_hat).real
        omega_y_phys = ifft2(omega_y_hat).real

        # Calculate nonlinear term (advection): u * grad(omega)
        advection_phys = u_phys * omega_x_phys + v_phys * omega_y_phys
        advection_hat = fft2(advection_phys)
        
        # Apply de-aliasing
        advection_hat *= dealias_mask

        # Diffusion term: nu * laplacian(omega)
        diffusion_hat = -nu * K2 * w_hat

        return diffusion_hat - advection_hat

    # --- 5. Time Stepping (Runge-Kutta 4) ---
    t = 0.0
    n_steps = int(end_time / dt)
    print_interval = int(0.5 / dt) # Print every 0.5 seconds

    print("Time      | Kinetic Energy | Enstrophy")
    print("-----------------------------------------")

    for i in range(n_steps + 1):
        # Calculate and print diagnostics
        if i % print_interval == 0:
            # Energy = 0.5 * sum(|u_hat|^2) = 0.5 * sum(|k|^{-2} |omega_hat|^2)
            energy = 0.5 * np.sum(K2_inv * np.abs(omega_hat)**2) / N**4
            # Enstrophy = 0.5 * sum(|omega_hat|^2)
            enstrophy = 0.5 * np.sum(np.abs(omega_hat)**2) / N**4
            print(f"{t:9.3f} | {energy:14.6e} | {enstrophy:12.6e}")

        if i == n_steps:
            break

        # RK4 steps
        k1 = compute_rhs(omega_hat)
        k2 = compute_rhs(omega_hat + 0.5 * dt * k1)
        k3 = compute_rhs(omega_hat + 0.5 * dt * k2)
        k4 = compute_rhs(omega_hat + dt * k3)
        
        omega_hat += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t += dt

if __name__ == '__main__':
    solve_and_observe_2d_navier_stokes()