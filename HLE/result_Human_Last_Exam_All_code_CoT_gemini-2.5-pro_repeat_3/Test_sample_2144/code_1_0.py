import numpy as np
from scipy.integrate import solve_ivp

def solve_trajectory():
    """
    Solves the given differential equation to find the position x0
    where the y-coordinate reaches -3.
    """
    
    # Define the differential equation dy/dx = f(x, y).
    # The equation is (dy/dx)^3 + y^2 = xy(dy/dx).
    # Let p = dy/dx. We need to solve the cubic p^3 - x*y*p + y^2 = 0 for p.
    def ode_func(x, y):
        # solve_ivp passes y as a single-element array.
        y_val = y[0]
        
        # Coefficients of the cubic polynomial in p: p^3 + 0*p^2 - xy*p + y^2 = 0
        coeffs = [1, 0, -x * y_val, y_val**2]
        
        # Find the roots of the cubic equation.
        roots = np.roots(coeffs)
        
        # For the domain of interest (x>=0, y<0), the cubic has a unique real root.
        # We select this real root as the slope p.
        p = roots[np.isreal(roots)].real[0]

        return p

    # Define an event function to find when the trajectory reaches y = -3.
    # The integration will stop when this function's value is zero.
    def reach_y_minus_3(x, y):
        return y[0] + 3

    # This makes the event terminal, so the integration stops.
    reach_y_minus_3.terminal = True
    # The event should trigger as y decreases through -3.
    reach_y_minus_3.direction = -1

    # Set the initial conditions from the problem statement.
    x_initial = 0
    y_initial = [-1]

    # Define the integration span for x. This should be large enough
    # for the event to occur.
    x_span = [x_initial, 20]

    # Use solve_ivp to find the solution.
    sol = solve_ivp(
        fun=ode_func,
        t_span=x_span,
        y0=y_initial,
        events=reach_y_minus_3,
        dense_output=True # Needed to get a continuous solution
    )

    # After integration, extract the result from the event that was triggered.
    if sol.t_events[0].size > 0:
        # The x-coordinate of the event
        x0 = sol.t_events[0][0]
        
        # The corresponding y-coordinate from the solution object
        y_final = sol.sol(x0)[0]

        # As requested, output the numbers in the final equation y(x0) = -3.
        print("The particle reaches the specified vertical coordinate.")
        print("The final equation is:")
        print(f"y({x0:.4f}) = {y_final:.0f}")

        # Return the value of x0 for the final answer.
        return x0
    else:
        print("The particle did not reach y = -3 within the integration range.")
        return None

# Run the solver and get the result.
x0_solution = solve_trajectory()

# Present the final answer in the required format.
if x0_solution is not None:
    print(f"<<<{x0_solution:.4f}>>>")
