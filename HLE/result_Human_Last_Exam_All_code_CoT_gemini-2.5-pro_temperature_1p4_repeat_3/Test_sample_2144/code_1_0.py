from scipy.integrate import solve_ivp
import numpy as np

def solve_trajectory():
    """
    Solves the particle's trajectory problem to find the position x0
    where y(x0) = -3.
    """
    
    # Define the ODE function dy/dx = f(x, y).
    # The function needs to solve for p = dy/dx from the cubic equation:
    # p^3 - x*y*p + y^2 = 0
    def ode_func(x, y):
        y_val = y[0]
        # Coefficients for the cubic polynomial in p
        coeffs = [1, 0, -x * y_val, y_val**2]
        # Find the roots of the cubic equation
        roots = np.roots(coeffs)
        
        # Along the trajectory, there is only one real root for p.
        # We select this real root as the value of dy/dx.
        real_root = roots[np.isreal(roots)].real[0]
        return real_root

    # Define an event function to stop the integration when y reaches -3.
    # The event occurs when this function returns zero.
    def event_y_equals_minus_3(x, y):
        return y[0] + 3
    
    # Make the event terminal, so the integration stops when it occurs.
    event_y_equals_minus_3.terminal = True
    # The direction of crossing is from positive (e.g., y=-1 -> y+3=2) to negative (y=-3.1 -> y+3=-0.1).
    event_y_equals_minus_3.direction = -1

    # Set the initial conditions
    x_initial = 0
    y_initial = [-1]

    # Define the integration span. We choose a range large enough for the event to occur.
    x_span = [0, 5]

    # Use solve_ivp to find the solution
    sol = solve_ivp(
        fun=ode_func,
        t_span=x_span,
        y0=y_initial,
        events=event_y_equals_minus_3,
        dense_output=True,
        rtol=1e-8, # Set a tighter tolerance for better accuracy
        atol=1e-8
    )

    # Extract the results from the solver output
    if sol.t_events[0].size > 0:
        x0 = sol.t_events[0][0]
        y0 = sol.y_events[0][0][0]
        
        print(f"The particle reaches the vertical coordinate y = -3 at position x0.")
        
        # As requested, output the numbers from the final equation y(x0) = -3
        print("\nThe final equation is y(x0) = -3. The numbers are:")
        print(f"x0 = {x0:.4f}")
        print(f"y = {y0:.0f}") # Should be exactly -3

    else:
        print("The particle did not reach y = -3 in the given integration interval.")

if __name__ == '__main__':
    solve_trajectory()