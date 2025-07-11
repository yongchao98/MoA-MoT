import numpy as np
from scipy.integrate import solve_ivp

def solve_trajectory():
    """
    Solves the de Broglie-Bohm trajectory problem.
    """
    # Step 1: Calculate the initial position x(0)
    # x(0) = 3 + (6(3-sqrt(3)))^(1/3) + (6(3+sqrt(3)))^(1/3)
    # This simplifies to 3 + (18 - 6*sqrt(3))^(1/3) + (18 + 6*sqrt(3))^(1/3)
    cbrt_term1 = (18 + 6 * np.sqrt(3))**(1/3)
    cbrt_term2 = (18 - 6 * np.sqrt(3))**(1/3)
    x0 = 3 + cbrt_term1 + cbrt_term2

    # Step 2: Define the ODE for the trajectory in the transformed frame x'(t)
    # dx'/dt = v'(x', t) = (4*x'*t) / ((x'^2 + t^2 - 3)^2 + 4*t^2)
    def velocity_transformed(t, x_prime):
        xp = x_prime[0]
        if t == 0:
            return [0.0]
        
        numerator = 4 * xp * t
        denominator = (xp**2 + t**2 - 3)**2 + 4 * t**2
        
        # This case should not be reached for the given trajectory
        if denominator == 0:
            return [np.inf]
            
        return [numerator / denominator]

    # Step 3: Set up and solve the ODE numerically
    xp0 = [x0]
    t_final = 2 * np.sqrt(3)
    t_span = [0, t_final]
    
    # Using a high-precision solver
    sol = solve_ivp(velocity_transformed, t_span, xp0, method='RK45', dense_output=True, rtol=1e-12, atol=1e-12)
    
    xp_final = sol.sol(t_final)[0]

    # Step 4: Transform the final position back to the original frame
    # x(t) = x'(t) - t^2 / 4
    x_final = xp_final - t_final**2 / 4

    # Print the equation and the final result
    print(f"The initial position is x(0) = {x0:.8f}")
    print(f"The target time is t = 2*sqrt(3) ≈ {t_final:.8f}")
    print(f"The position in the transformed frame at the target time is x'(t) ≈ {xp_final:.8f}")
    print(f"The final position x(t) is calculated by x(t) = x'(t) - t^2 / 4")
    print(f"So, x({t_final:.8f}) = {xp_final:.8f} - ({t_final:.8f})^2 / 4")
    print(f"x({t_final:.8f}) = {xp_final:.8f} - {t_final**2:.8f} / 4")
    print(f"x({t_final:.8f}) = {xp_final:.8f} - {t_final**2 / 4:.8f}")
    print(f"The final result is x({t_final:.8f}) ≈ {x_final:.8f}")

solve_trajectory()