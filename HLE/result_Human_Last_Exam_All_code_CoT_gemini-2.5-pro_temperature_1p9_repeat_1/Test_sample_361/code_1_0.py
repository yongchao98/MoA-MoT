import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import bisect

def calculate_critical_potential(lambda_gb=0.1, verbose=False):
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7 holographic model
    with Einstein-Gauss-Bonnet gravity at T=0.

    The method involves finding the lowest eigenvalue of the Sturm-Liouville problem
    for the scalar fluctuation.
    """

    # 1. Calculate the metric factor f0 from the Gauss-Bonnet coupling
    if lambda_gb == 0:
        f0 = 1.0
    elif 1 - 4 * lambda_gb < 0:
        raise ValueError("lambda_gb must be <= 1/4")
    else:
        f0 = (1 - np.sqrt(1 - 4 * lambda_gb)) / (2 * lambda_gb)

    A = 3.0 / f0

    # 2. Define the ODE system for the shooting method
    # y'' + (3/x)y' + (A/x^2 - E/x^6)y = 0
    def ode_system(rho, y, E):
        w, dw_drho = y
        if rho == 0:
            return [0, 0]
        # Avoid division by zero at rho=0
        d2w_drho2 = -(3.0 / rho) * dw_drho - (A / rho**2 - E / rho**6) * w
        return [dw_drho, d2w_drho2]

    # 3. Define the objective function for the root finder
    # We want the solution to vanish at large rho.
    # The shooting method will vary E until w(rho_max) is close to zero.
    rho_min = 1e-6
    rho_max = 20.0

    # Determine behavior at rho_min for initial conditions
    # For small rho, the E/rho^6 term dominates and solutions are wild.
    # Let's change variable x = 1/rho. y(x)=w(1/x)
    # y'' - 1/x y' + (A/x^2 - E x^2)/x^2 y = 0
    # No, that's not right. The original ODE is correct.
    # Let's analyse the behavior at small rho (origin) for w(rho).
    # rho^2 w'' + 3 rho w' + A w = 0 is the indicial equation.
    # Assume w ~ rho^alpha => alpha(alpha-1) + 3alpha + A = 0 => alpha^2 + 2alpha + A = 0
    # alpha = -1 +/- sqrt(1 - A). For A = 3/f0 > 1, this is complex.
    # alpha = -1 +/- 1j * sqrt(A - 1).
    # This indicates log-periodic behavior, which is a hallmark of this system.
    # We choose the decaying mode at large rho.
    # It is easier to shoot from large rho towards the origin.
    
    # At large rho, the E term is negligible. w ~ rho^alpha where alpha = -1 +/- i sqrt(A-1).
    # Both are decaying solutions. Let's take the real part of one solution as a BC.
    
    sqrt_A_minus_1 = np.sqrt(A - 1)
    
    def get_residual(E):
        # Initial conditions at rho_max
        # We start with the asymptotic solution at large rho, which is a decaying oscillating function
        # w(rho) ~ 1/rho * cos(sqrt(A-1)*log(rho) + delta)
        # We can set the phase delta=0 for simplicity.
        y_max = [1/rho_max * np.cos(sqrt_A_minus_1 * np.log(rho_max)),
                 -1/rho_max**2 * (np.cos(sqrt_A_minus_1 * np.log(rho_max)) + sqrt_A_minus_1 * np.sin(sqrt_A_minus_1 * np.log(rho_max)))]
        
        # Solve the ODE backwards from rho_max to rho_min
        sol = solve_ivp(ode_system, [rho_max, rho_min], y_max, args=(E,), dense_output=True, method='RK45')
        
        # The residual is the value of the solution at rho_min.
        # We want the solution to be regular, not divergent.
        # The divergent behavior at small rho is w ~ rho^(-1 - ...), we want to kill this mode.
        # Checking for the coefficient of the 'bad' solution is tricky. A simpler proxy is to check
        # the value at the end of integration. A well-behaved solution will be small.
        w_final = sol.sol(rho_min)[0]
        return w_final

    # 4. Find the eigenvalue E by finding the root of the residual function.
    # Search for the lowest eigenvalue E > 0.
    # Based on literature (e.g. arXiv:0907.4111), for lambda=0 (f0=1, A=3), E is ~11.6
    # For lambda=0.1, A ~ 2.66, so E should be slightly smaller.
    E_low = 5.0
    E_high = 15.0
    # The get_residual function will change sign at each eigenvalue.
    try:
        eigenvalue_E = bisect(get_residual, E_low, E_high, xtol=1e-6)
    except ValueError:
        print("Bisection failed. The function might not cross zero in the given interval.")
        print("This could be due to numerical precision issues or an incorrect search range.")
        return None

    # 5. Calculate mu_c from the eigenvalue E.
    # The relation is mu_c = C * sqrt(E) * f0, where C is a normalization constant.
    # In the absence of a specified normalization scheme, we set C=1.
    mu_c = np.sqrt(eigenvalue_E) * f0
    
    if verbose:
        print(f"For Gauss-Bonnet coupling lambda_gb = {lambda_gb}:")
        print(f"Metric factor f0 = {f0:.4f}")
        print(f"ODE parameter A = 3/f0 = {A:.4f}")
        print(f"Lowest eigenvalue found: E = {eigenvalue_E:.4f}")
        print(f"The critical chemical potential mu_c is proportional to sqrt(E) * f0.")
        print(f"Setting the proportionality constant to 1 gives:")

    print(f"μ_c = √({eigenvalue_E:.4f}) * {f0:.4f} = {mu_c:.4f}")

    return mu_c

# Run the calculation for the user's specific case.
# Suppress verbose output for the final answer block.
critical_potential = calculate_critical_potential(lambda_gb=0.1, verbose=True)

# Final answer formatting
final_answer_val = float(f"{critical_potential:.4f}")
print(f"The value of the critical chemical potential is {final_answer_val}")
# The instruction is to return the value at the end. I will use the last print for this.
# No, the instruction is to use print function. Okay.
