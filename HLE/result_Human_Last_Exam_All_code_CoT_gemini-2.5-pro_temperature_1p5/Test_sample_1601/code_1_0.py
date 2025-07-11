import numpy as np
from scipy.integrate import solve_ivp
import warnings

# Suppress warnings that may arise from overflows during integration,
# as these are expected signs of solution blow-up.
warnings.filterwarnings("ignore", category=RuntimeWarning)

def solve_differential_equation_and_estimate_area():
    """
    Solves the given system of ODEs for a grid of initial conditions
    and estimates the area of the set Omega that leads to a specific blow-up behavior.
    """

    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -b**2 / 2 - np.exp(t) * a**2 - a
        return [dadt, dbdt]

    # Define the domain for the initial conditions [a(0), b(0)]
    a0_min, a0_max = -10, 1
    b0_min, b0_max = 10, 20

    # Create a grid of initial points to sample the domain.
    # A finer grid provides a more accurate estimate.
    a_step, b_step = 0.25, 0.25
    a0_vals = np.arange(a0_min, a0_max + a_step, a_step)
    b0_vals = np.arange(b0_min, b0_max + b_step, b_step)

    total_points = len(a0_vals) * len(b0_vals)
    blowup_count = 0

    # Set the time span for the numerical integration.
    # Due to the exp(t) term, any blow-up is expected to occur relatively quickly.
    t_span = [0, 5]
    t_max = t_span[1]

    # Iterate through each initial condition on the grid
    for a0 in a0_vals:
        for b0 in b0_vals:
            y0 = [a0, b0]
            
            # Solve the ODE system for the current initial condition.
            sol = solve_ivp(ode_system, t_span, y0, method='RK45', atol=1e-6, rtol=1e-3)

            # Check if the solution blew up as expected.
            # This is primarily indicated by the solver failing to reach t_max.
            if sol.t[-1] < t_max:
                a_final, b_final = sol.y[0, -1], sol.y[1, -1]
                
                # The specified blow-up (a->inf, b->-inf) happens if the trajectory
                # enters the a>0, b<0 quadrant.
                if a_final > 0 and b_final < 0:
                    blowup_count += 1
    
    # Calculate the total area of the initial conditions domain
    total_area = (a0_max - a0_min) * (b0_max - b0_min)
    
    # Estimate the measure of the set Omega
    if total_points > 0:
        m_omega_estimate = (blowup_count / total_points) * total_area
    else:
        m_omega_estimate = 0

    print("Numerical Estimation Results:")
    print(f"Total points sampled in domain D: {total_points}")
    print(f"Points in Omega (leading to blow-up): {blowup_count}")
    print(f"Total area of domain D: {total_area}")
    print(f"Estimated area of Omega, m(Omega): {m_omega_estimate:.4f}")
    
    # Outputting the components of the final calculation as requested
    print("\nFinal area estimation equation:")
    print(f"m(Omega) = ({blowup_count} / {total_points}) * {total_area}")


solve_differential_equation_and_estimate_area()