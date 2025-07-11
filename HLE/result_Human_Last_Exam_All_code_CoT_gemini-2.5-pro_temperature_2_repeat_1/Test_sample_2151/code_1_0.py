import numpy as np
from scipy.integrate import solve_ivp
import warnings

# Suppress potential RuntimeWarning that can occur with np.exp on longdouble types
warnings.filterwarnings("ignore", category=RuntimeWarning)

def solve_pde():
    """
    This function solves the given PDE numerically and prints the result.
    """
    # Step 1: Define simulation parameters
    L = 40.0       # The spatial domain will be [-L, L]
    N = 4001       # Number of grid points (odd to include x=0)
    x = np.linspace(-L, L, N, dtype=np.longdouble)
    dx = x[1] - x[0]
    x_0_index = N // 2  # The index of the grid point for x=0

    # Step 2: Define and evaluate the initial conditions on the grid
    # u(x, 0)
    def u_initial(x_grid):
        ex = np.exp(x_grid)
        return -2.0 + (1.0 - np.tanh(x_grid)) / (ex + 1.0)

    # ∂u/∂t(x, 0)
    def ut_initial(x_grid):
        tanh_x = np.tanh(x_grid)
        sech_x = 1.0 / np.cosh(x_grid)
        sech_half_x_sq = (1.0 / np.cosh(x_grid / 2.0))**2
        return 0.25 * (tanh_x - 1.0) * sech_half_x_sq * (tanh_x - sech_x - 2.0)

    # Set up the initial state vector y = [u_0, ..., u_{N-1}, v_0, ..., v_{N-1}] where v=du/dt
    u0 = u_initial(x).astype(np.float64)
    v0 = ut_initial(x).astype(np.float64)
    y0 = np.concatenate((u0, v0))

    # Step 3: Define the ODE system for the Method of Lines
    def pde_system(t, y):
        u = y[:N]
        v = y[N:]
        
        # Initialize time derivatives. Boundary values are fixed, so their derivatives are 0.
        dudt = np.zeros_like(u)
        dvdt = np.zeros_like(v)
        
        # For interior points, du/dt = v
        dudt[1:N-1] = v[1:N-1]
        
        # Calculate spatial derivatives for interior points using central differences
        u_int = u[1:N-1]
        du_dx = (u[2:] - u[:-2]) / (2 * dx)
        d2u_dx2 = (u[2:] - 2*u[1:-1] + u[:-2]) / (dx**2)
        
        # The PDE is (1/8)u_tt + u_t + F(u) = 0, where F(u) = u*u_x - (1/8)u_xx - (u-1)u(u+2)
        # Let v=u_t, then dv/dt = u_tt. The ODE for v is dv/dt = -8*(v + F(u))
        F_u = u_int * du_dx - (1.0/8.0) * d2u_dx2 - (u_int - 1.0) * u_int * (u_int + 2.0)
        dvdt[1:N-1] = -8.0 * (v[1:N-1] + F_u)
        
        return np.concatenate((dudt, dvdt))

    # Step 4: Solve the system of ODEs
    print("Starting numerical integration. This may take a moment...")
    t_span = [0.0, 1.0]
    # We use a solver suitable for stiff problems, like 'Radau'.
    sol = solve_ivp(pde_system, t_span, y0, method='Radau', t_eval=[1.0])
    print("Numerical integration complete.")

    # Step 5: Extract the result and print the final calculation
    if sol.status == 0:
        y_final = sol.y[:, -1]
        u_final = y_final[:N]
        u_0_1 = u_final[x_0_index]
        result = -u_0_1 / 2.0
        
        print("\nFinal Calculation:")
        print(f"The value of u at x=0 and t=1 is found to be u(0,1) = {u_0_1}")
        print("The quantity to determine is -u(0,1)/2. The final equation with the computed values is:")
        print(f"-({u_0_1}) / 2 = {result}")
    else:
        print("The ODE solver failed.")
        print(sol.message)

if __name__ == '__main__':
    solve_pde()