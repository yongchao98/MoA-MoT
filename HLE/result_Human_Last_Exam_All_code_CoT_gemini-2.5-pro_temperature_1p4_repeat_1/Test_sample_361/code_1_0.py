import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def calculate_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using the shooting method.
    """
    # 1. Define the model parameters
    lambda_gb = 0.1  # Given Gauss-Bonnet coupling
    m_squared = -3.0 # Mass squared of the scalar field for the D3/D7 setup (m^2 * L^2)
    q_charge = 1.0   # Charge of the scalar field (can be set to 1)

    # Set the scale by fixing the horizon radius. The temperature is T = r_h / pi.
    r_h = 1.0

    # 2. Define the EGB metric function f(r) and its derivative
    # f(r) = (r^2 / (2*lambda)) * [1 - sqrt(1 - 4*lambda*(1 - r_h^4/r^4))]
    def f(r):
        r_safe = np.maximum(r, r_h + 1e-9) # Avoid issues at r=r_h
        term_under_sqrt = 1 - 4 * lambda_gb * (1 - (r_h**4 / r_safe**4))
        return (r**2 / (2 * lambda_gb)) * (1 - np.sqrt(term_under_sqrt))

    def f_prime(r):
        if r <= r_h: # Handle the limit at the horizon
            return 4 * r_h
        r_safe = np.maximum(r, r_h + 1e-9)
        g = np.sqrt(1 - 4 * lambda_gb + 4 * lambda_gb * (r_h**4 / r_safe**4))
        return r / lambda_gb * (1 - g) + 4 * r_h**4 / (r_safe**3 * g)

    # 3. Define the function for the ODE system to be solved
    # System: y[0] = psi, y[1] = psi'
    def ode_system(r, y, mu):
        psi, psi_prime = y
        
        fr = f(r)
        if fr < 1e-9: # Should not be called at the horizon itself
            return [0, 0]
        fpr = f_prime(r)
        
        # A_t(r) = mu * (1 - r_h^2 / r^2)
        At = mu * (1 - r_h**2 / r**2)
        
        # The ODE: psi'' = - (3/r + f'/f) * psi' - (m^2/f + q^2*A_t^2/f^2) * psi
        d_psi_prime = -(3/r + fpr/fr) * psi_prime - (m_squared/fr + (q_charge**2 * At**2)/fr**2) * psi
        return [psi_prime, d_psi_prime]

    # 4. Define the target function for the root finder
    # This function solves the ODE for a given mu and checks the boundary condition at infinity.
    def solve_and_check_bc(mu, r_start, r_end):
        # Initial conditions determined by regularity at the horizon r_h
        psi_start = 1.0
        psi_prime_start = -m_squared / (4 * r_h**3)
        y0 = [psi_start, psi_prime_start]
        
        # Solve the ODE
        sol = solve_ivp(
            fun=lambda r, y: ode_system(r, y, mu),
            t_span=[r_start, r_end],
            y0=y0,
            method='RK45'
        )
        
        # We need the solution to be normalizable, i.e., have the correct fall-off at infinity.
        # This can be checked by evaluating a specific combination of psi and psi' at r_end -> infinity
        # that should be zero for the correct solution.
        alpha = (1 - np.sqrt(1 - 4*lambda_gb))/(2*lambda_gb)
        delta_plus = 2 + np.sqrt(4 + m_squared/alpha)
        
        psi_end = sol.y[0][-1]
        psi_prime_end = sol.y[1][-1]
        
        # The target function that is zero for the correct mu
        target_value = r_end * psi_prime_end + delta_plus * psi_end
        return target_value

    # 5. Execute the shooting method
    r_start_integration = r_h + 1e-5
    r_end_integration = 200.0  # A large number to approximate infinity
    
    # Bracket for the chemical potential. Based on existing literature, the value is of order 1.
    mu_bracket = [1.0, 10.0]

    try:
        # Use a robust root-finding algorithm
        result = root_scalar(
            lambda mu: solve_and_check_bc(mu, r_start_integration, r_end_integration),
            bracket=mu_bracket,
            method='brentq'
        )
        critical_mu = result.root

        # Print the parameters used and the final result
        print("Calculation of the critical chemical potential (μ_c)")
        print("-" * 50)
        print(f"Holographic model: 5D Einstein-Gauss-Bonnet on D3/D7 background")
        print(f"Parameters used in the final calculation:")
        print(f"  Gauss-Bonnet coupling (λ_GB) = {lambda_gb}")
        print(f"  Scalar field mass squared (m²) = {m_squared}")
        print(f"  Scalar field charge (q)       = {q_charge}")
        print("-" * 50)
        print(f"The resulting critical chemical potential is μ_c = {critical_mu:.4f}")
        print(f"\n<<<__{critical_mu:.4f}__>>>")

    except ValueError:
        print("Could not find the critical potential. The bracket might be incorrect.")

if __name__ == '__main__':
    calculate_critical_potential()