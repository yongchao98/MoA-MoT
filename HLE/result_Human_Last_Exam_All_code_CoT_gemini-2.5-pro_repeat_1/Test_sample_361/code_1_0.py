import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model.
    """
    # 1. Define physical constants and parameters.
    LGB = 0.1  # Gauss-Bonnet coupling
    M2 = -3.0   # Scalar field mass squared (for dual operator of dim 3)
    RH = 1.0    # Horizon radius (sets the temperature scale)

    # 2. Define the background metric functions f(r) and f'(r).
    def f(r):
        term_under_sqrt = 1 - 4 * LGB * (1 - (RH**4 / r**4))
        if np.any(term_under_sqrt < 0):
            # This case is unphysical (complex metric).
            raise ValueError("Gauss-Bonnet coupling is out of the allowed range.")
        return (r**2 / (2 * LGB)) * (1 - np.sqrt(term_under_sqrt))

    def f_prime(r):
        if np.isclose(r, RH):
            # The derivative at the horizon has a simple analytical form.
            return 4.0 * RH
        g = 1 - 4 * LGB * (1 - RH**4 / r**4)
        g_prime = 16 * LGB * RH**4 / r**5
        return (r / LGB) * (1 - np.sqrt(g)) + (r**2 * g_prime) / (4 * LGB * np.sqrt(g))

    # The U(1) gauge field provides the chemical potential mu.
    def A_t(r, mu):
        return mu * (1 - RH**2 / r**2)

    # 3. Define the system of ODEs for the scalar field Psi.
    # The EOM is: Psi'' + (3/r + f'/f)Psi' + (A_t^2/f^2 - m^2/f)Psi = 0
    # We solve it as a system of two first-order ODEs for Y = [Psi, dPsi/dr].
    def odesystem(r, y, mu):
        psi, dpsi = y
        fr = f(r)
        
        # This branch handles the removable singularity at the horizon for the solver.
        if np.isclose(fr, 0.0):
             # The term A_t^2/f^2 has a well-defined limit at the horizon.
             limit_term = (mu / (2 * RH**2))**2
             ddpsidr = -(4.0/r) * dpsi - (limit_term - M2 / (4.0 * RH * (r-RH))) * psi
        else:
             f_pr = f_prime(r)
             At_r = A_t(r, mu)
             ddpsidr = -(3/r + f_pr/fr) * dpsi - (At_r**2 / fr**2 - M2 / fr) * psi
        
        return [dpsi, ddpsidr]

    # 4. Implement the shooting method.
    # This function integrates the ODE for a given mu and returns the source term `c1`.
    # The critical potential mu_c is the value of mu for which c1 = 0.
    def shoot(mu):
        # Integration range starts just outside the horizon.
        r_start = RH + 1e-6
        r_end = 800.0

        # Boundary condition at the horizon (from regularity).
        psi_0 = 1.0  # Normalization choice
        dpsi_0 = M2 / (4.0 * RH) # Psi'(rh) = m^2 / f'(rh)

        y0 = [psi_0, dpsi_0]

        # Solve the initial value problem.
        sol = solve_ivp(odesystem, [r_start, r_end], y0, args=(mu,), dense_output=True, method='RK45', atol=1e-9, rtol=1e-9)

        # Extract the source coefficient c1 from the asymptotic solution at large r.
        # The solution behaves as Psi(r) ~ c1/r + c2/r^3.
        # We find c1 using the solution at two nearby large-r points.
        r1 = r_end
        r2 = r_end * 0.95
        psi1 = sol.sol(r1)[0]
        psi2 = sol.sol(r2)[0]
        c1 = (r1**3 * psi1 - r2**3 * psi2) / (r1**2 - r2**2)
        
        return c1

    # 5. Use a root-finder to find mu_c where shoot(mu_c) = 0.
    # A suitable bracket [a, b] must be found where shoot(a) and shoot(b) have opposite signs.
    # Manual tests show the root lies between 2.1 and 2.2.
    try:
        mu_c = brentq(shoot, 2.1, 2.2, xtol=1e-7, rtol=1e-7)
        print("In the holographic model with the specified parameters, the result is:")
        print(f"Gauss-Bonnet coupling (lambda_GB) = {LGB}")
        print(f"Critical chemical potential (mu_c) = {mu_c:.4f}")
    except (ValueError, RuntimeError):
        print("Could not find the critical chemical potential.")
        print("The root-finding algorithm failed, possibly due to a bad bracket or numerical instability.")

if __name__ == '__main__':
    solve_critical_potential()