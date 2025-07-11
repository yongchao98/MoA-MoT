import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

# Plan:
# 1. Define the parameters of the holographic model.
#    - Gauss-Bonnet coupling lambda_gb = 0.1
#    - Operator dimension Delta = 3, which corresponds to scalar mass squared m^2 = Delta * (Delta - d) = 3 * (3-4) = -3. (Bulk dimension is D=5, so boundary is d=4).
# 2. Set up the differential equation for the scalar field perturbation psi at the critical point.
#    - The equation is a second-order linear ODE for psi(z), where z is the holographic coordinate (z=0 is boundary, z=1 is horizon).
#    - The equation involves a function f(z) from the EGB black hole metric and the chemical potential mu, which will be our eigenvalue.
# 3. Formulate the problem as a shooting method problem.
#    - We need to find the value of mu (mu_c) for which a normalizable solution for psi(z) exists.
#    - We integrate the ODE from the horizon (z=1) to the boundary (z=0).
#    - At the horizon (z=1), we impose regularity, which sets the initial conditions for psi and its derivative psi'.
#    - At the boundary (z=0), the normalizability condition (source-free) requires that the solution behaves as z^3. This gives us a condition to check.
# 4. Numerically solve for mu_c.
#    - Use `solve_ivp` to integrate the ODE for a guessed mu.
#    - Use a root-finding algorithm (`root_scalar`) to find the value of mu that satisfies the boundary condition at z=0.
# 5. Print the final result and the equation solved.

# --- Model Parameters ---
lambda_gb = 0.1
m2 = -3  # Mass squared for a scalar dual to an operator with dimension Delta=3 in d=4

# --- Metric Function f(z) and its derivative (numerically stable forms) ---
# f(z) for the 5D EGB-AdS black hole metric. z = r_h / r.
# We set AdS radius L=1 and horizon radius r_h=1.
def f_stable(z, lam):
    if z == 1.0:
      return 0.0
    # This form avoids sqrt(small_number) for z~1, improving stability
    return (2*(1-z**4))/(1+np.sqrt(1-4*lam*(1-z**4)))
    
def dfdz_stable(z,lam):
    if z==0.0:
        return 0.0
    sq = np.sqrt(1-4*lam*(1-z**4))
    if z == 1.0: # l'Hopital's rule on f'(z) shows it's -4
        return -4.0
    term1 = -8*z**3 / (1+sq)
    term2 = (2*(1-z**4) * (8*lam*z**3) / sq) / (1+sq)**2
    return term1 - term2

# --- ODE System ---
# Equation for psi: psi'' + p(z)psi' + q(z)psi = 0
# p(z) = f'(z)/f(z) + 3/z
# q(z) = mu^2 * (1-z^2)^2 / f(z)^2 - m^2 / (z^2 * f(z))
def odesystem(z, y, mu, lam, m2_val):
    psi, dpsi = y
    
    # Use stable versions of f and dfdz
    fz = f_stable(z, lam)
    dfz = dfdz_stable(z, lam)
    
    p_z = (dfz / fz) + (3 / z)
    
    term_mu = (mu**2 * (1 - z**2)**2) / (fz**2)
    term_m2 = -m2_val / (z**2 * fz)
    
    q_z = term_mu + term_m2
    
    # y[0] = psi, y[1] = dpsi/dz
    # d(psi)/dz = dpsi
    # d(dpsi)/dz = -p(z)*dpsi - q(z)*psi
    return [dpsi, -p_z * dpsi - q_z * psi]

# --- Shooting Method ---
def objective_function(mu, lam, m2_val):
    """
    Solves the ODE and returns a value that should be zero when mu is the correct eigenvalue.
    """
    # Integration range, away from singular points z=0 and z=1
    z_start = 1.0 - 1e-6
    z_end = 1e-6
    
    # Initial conditions at the horizon (z=1).
    # Regularity requires psi'(1) = 0 for our EOM, and we can set psi(1)=1 by scaling.
    y_start = [1.0, 0.0]
    
    # Solve the IVP
    sol = solve_ivp(
        odesystem, 
        [z_start, z_end], 
        y_start, 
        args=(mu, lam, m2_val),
        dense_output=True,
        method='RK45'
    )
    
    # Extract solution at the end point z_end
    psi_end, dpsi_end = sol.sol(z_end)
    
    # Boundary condition at z=0. The solution should be of the form C * z^3.
    # The other solution is z^1.
    # To enforce the coefficient of z^1 to be zero, we can check a Wronskian-like quantity.
    # If psi(z) = A*z^1 + B*z^3, then 3*psi(z) - z*psi'(z) = 2*A*z. We want A=0.
    target_value = 3 * psi_end - z_end * dpsi_end
    
    return target_value

# --- Main Execution ---
print("Attempting to find the critical chemical potential μ_c...")
# Bracket for the root finder. Based on literature, the value should be around 5.
bracket = [4.0, 6.0]

try:
    # Find the root mu_c
    result = root_scalar(
        objective_function, 
        args=(lambda_gb, m2), 
        bracket=bracket,
        method='brentq'
    )
    mu_c = result.root

    print("\nCalculation successful.\n")

    print("The differential equation being solved for the scalar field ψ(z) is:")
    print("ψ''(z) + [f'(z)/f(z) + 3/z] ψ'(z) + [μ²(1-z²)²/f(z)² - m²/ (z²f(z))] ψ(z) = 0\n")

    print(f"For the given parameters:")
    print(f"Gauss-Bonnet coupling, λ_GB = {lambda_gb}")
    print(f"Operator dimension, Δ = 3, which implies mass squared, m² = {m2}\n")

    print("The numerical calculation gives:")
    print(f"Critical chemical potential, μ_c = {mu_c:.4f}")
    
    T = 1 / np.pi
    mu_c_over_T = mu_c * np.pi
    print(f"(For a black hole temperature T = 1/π ≈ {T:.4f})")
    print(f"The dimensionless ratio, μ_c/T ≈ {mu_c_over_T:.4f}\n")
    
    # Print the numbers in the final equation as requested for a sample z=0.5
    z_sample = 0.5
    f_sample = f_stable(z_sample, lambda_gb)
    df_sample = dfdz_stable(z_sample, lambda_gb)
    coeff_psi_prime = (df_sample / f_sample) + (3 / z_sample)
    coeff_psi_mu_part = (1 - z_sample**2)**2 / f_sample**2
    coeff_psi_m2_part = -m2 / (z_sample**2 * f_sample)
    coeff_psi = coeff_psi_mu_part * mu_c**2 + coeff_psi_m2_part

    print("Example coefficients of the equation at z=0.5 with the calculated μ_c:")
    # final equation is psi'' + coeff_prime * psi' + coeff_psi * psi = 0
    # For printing, we show the numbers that make up the final coefficient
    print(f"1.0000 * ψ''(0.5) + {coeff_psi_prime:.4f} * ψ'(0.5) + ({coeff_psi_mu_part:.4f}*{mu_c:.4f}² + {coeff_psi_m2_part:.4f}) * ψ(0.5) = 0")
    print(f"which evaluates to:")
    print(f"1.0000 * ψ''(0.5) + {coeff_psi_prime:.4f} * ψ'(0.5) + {coeff_psi:.4f} * ψ(0.5) = 0")
    
except (ValueError, RuntimeError) as e:
    print(f"Error during calculation: {e}")
    print("The bracket might not contain a root, or there was a numerical issue.")
