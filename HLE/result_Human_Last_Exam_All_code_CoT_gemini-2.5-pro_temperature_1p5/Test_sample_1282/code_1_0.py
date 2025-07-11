import numpy as np
from scipy.integrate import solve_ivp
import warnings

# Suppress ComplexWarning that may arise from numerical instabilities
warnings.filterwarnings("ignore", category=np.ComplexWarning)

def demonstrate_blowup():
    """
    Demonstrates the blow-up behavior of the equation ∂t u + (1+t)Δu = 0
    by solving a 1D version: ∂t u + (1+t)∂xx u = 0.
    """
    # 1. Define the spatial domain and discretization
    L = 2 * np.pi  # Length of the domain
    N = 128        # Number of spatial points
    x = np.linspace(0, L, N, endpoint=False) # Grid points
    dx = x[1] - x[0]

    # 2. Define the initial condition (a smooth function)
    k = 2 # Wavenumber
    u0 = np.sin(k * x)

    # 3. Define the ODE system resulting from spatial discretization
    # We use a finite difference approximation for the second derivative (∂xx u)
    # with periodic boundary conditions. The resulting matrix is the discrete Laplacian.
    # The ODE is d(u_vec)/dt = -(1+t) * Laplacian * u_vec
    
    # Construct the Laplacian matrix D2
    main_diag = -2 * np.ones(N)
    off_diag = np.ones(N - 1)
    D2 = (np.diag(main_diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)) / dx**2
    # Add periodic boundary conditions
    D2[0, -1] = 1 / dx**2
    D2[-1, 0] = 1 / dx**2

    def backward_heat_equation(t, u):
        """The system of ODEs: du/dt = f(t, u)"""
        return -(1 + t) * (D2 @ u)

    # 4. Solve the ODE system
    # We solve over a very short time interval because the solution blows up quickly.
    t_span = [0, 0.4]
    t_eval = np.linspace(t_span[0], t_span[1], 10)

    solution = solve_ivp(
        backward_heat_equation,
        t_span,
        u0,
        t_eval=t_eval,
        method='RK45'
    )
    
    print("--- Numerical Demonstration of Blow-up ---")
    print(f"Simulating ∂t u + (1+t)∂xx u = 0 for initial data u(x,0)=sin({k}x).\n")
    print("The anti-dissipative term +(1+t)Δu causes the solution's norm to grow explosively.")
    print("We verify the energy evolution: d/dt(||u||^2) = 2(1+t) * ||∂x u||^2.\n")
    print(f"{'Time (t)':<10} | {'L2 Norm ||u||':<18} | {'d/dt(||u||^2)':<18} | {'2(1+t)||∂x u||^2':<18}")
    print("-" * 75)

    # 5. Analyze and print the results
    prev_norm_sq = None
    prev_t = None
    for i, t in enumerate(solution.t):
        u = solution.y[:, i]
        
        # Calculate terms of the energy evolution equation
        # Note: integration is approximated by summation * dx
        
        # ||u||^2 = ∫|u|^2 dx
        norm_sq = np.sum(u**2) * dx
        
        # ||∂x u||^2 = ∫|∂x u|^2 dx
        grad_u = np.gradient(u, dx)
        enstrophy = np.sum(grad_u**2) * dx
        
        # Left-hand side: d/dt(||u||^2)
        if i > 0:
            dt = t - prev_t
            d_norm_sq_dt = (norm_sq - prev_norm_sq) / dt
        else:
            d_norm_sq_dt = np.nan # Cannot compute derivative at t=0
        
        # Right-hand side: 2 * (1+t) * ||∂x u||^2
        rhs_energy_eq = 2 * (1 + t) * enstrophy

        print(f"{t:<10.4f} | {np.sqrt(norm_sq):<18.4e} | {d_norm_sq_dt:<18.4e} | {rhs_energy_eq:<18.4e}")

        prev_norm_sq = norm_sq
        prev_t = t
        # Stop if the norm becomes excessively large
        if np.sqrt(norm_sq) > 1e10:
            print("\nNorm exceeded 1e10. Stopping simulation as blow-up is evident.")
            break

if __name__ == '__main__':
    demonstrate_blowup()