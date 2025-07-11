import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma

def solve_balloon_radius():
    """
    Solves the given differential equation for the balloon radius y(t) at t=pi/4.
    """
    # Calculate the initial condition y(0)
    # y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1)
    y0_val = 1.0 / (128 * np.power(3, 1/6) * gamma(2/3))

    # Define the system of first-order ODEs from the given second-order ODE.
    # Y[0] = y(t), Y[1] = y'(t) = v(t)
    def ode_system(t, Y):
        y, v = Y
        
        sint = np.sin(t)
        
        # The ODE is singular at t=0, so we should not evaluate it there.
        # Also, the RHS has sqrt(sin(t)), so sin(t) must be non-negative.
        if t == 0 or sint < 0:
            return [0, 0] # Should not be reached if t_start > 0
            
        cost = np.cos(t)
        tant = sint / cost
        sect = 1 / cost
        
        # Coefficients of the ODE in the form: C2*y'' + C1*y' + C0*y = R_term
        t4_plus_1 = t**4 + 1
        
        C2 = 4 * t4_plus_1 * tant * sect
        C1 = 2 * (t**4 + 2 * tant * (t4_plus_1 * tant + 8 * t**3) + 1) * sect
        C0 = 8 * t**2 * (t + 2 * tant * (t * tant + 3)) * sect
        R_term = t4_plus_1 * y * np.sqrt(sint)
        
        # y'' = v' = (R_term - C1*v - C0*y) / C2
        dvdt = (R_term - C1 * v - C0 * y) / C2
        dydt = v
        
        return [dydt, dvdt]

    # The ODE is singular at t=0. We start integration at a small t_start > 0.
    t_start = 1e-8
    t_end = np.pi / 4

    # Determine initial conditions at t_start using asymptotic analysis near t=0.
    # y(t) ≈ y(0) + (y(0)/9) * t^(3/2)
    # y'(t) ≈ (y(0)/6) * t^(1/2)
    y_start = y0_val + (y0_val / 9) * np.power(t_start, 1.5)
    v_start = (y0_val / 6) * np.power(t_start, 0.5)

    initial_conditions = [y_start, v_start]
    time_span = [t_start, t_end]

    # Solve the ODE using a robust method (RK45) with tight tolerances.
    solution = solve_ivp(ode_system, time_span, initial_conditions, method='RK45', dense_output=True, rtol=1e-9, atol=1e-12)

    # Evaluate the solution at the final time t = pi/4
    y_pi_over_4 = solution.sol(t_end)[0]

    print(f"{y_pi_over_4}")

solve_balloon_radius()