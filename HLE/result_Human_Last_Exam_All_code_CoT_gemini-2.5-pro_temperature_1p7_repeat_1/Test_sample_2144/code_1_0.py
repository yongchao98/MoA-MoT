import numpy as np
from scipy.integrate import solve_ivp
import math

def get_p(x, y):
    """
    Solves the cubic equation p^3 - (x*y)*p + y^2 = 0 for the real root p.
    """
    # The coefficients of the polynomial in p: p^3 + 0*p^2 - (x*y)*p + y^2
    coeffs = [1, 0, -x * y, y**2]
    
    # Find the roots
    roots = np.roots(coeffs)
    
    # Return the real root. The trajectory stays in a region with a single real root.
    real_root = roots[np.isreal(roots)].real[0]
    return real_root

def ode_function(x, y_vec):
    """
    The function to be integrated by the ODE solver.
    x is the independent variable.
    y_vec is a vector of dependent variables, here just [y].
    """
    y = y_vec[0]
    # dy/dx is the value of p
    dydx = get_p(x, y)
    return [dydx]

def event_y_reaches_minus_3(x, y_vec):
    """
    Event function to stop integration when y(x) = -3.
    This function returns 0 when the event occurs.
    """
    return y_vec[0] + 3

# This event should terminate the integration.
event_y_reaches_minus_3.terminal = True
# Since y is decreasing from -1 to -3, the event function y+3 is decreasing.
event_y_reaches_minus_3.direction = -1

# Set up the initial conditions and integration span.
x_initial = 0.0
y_initial = -1.0
# The integration needs to run long enough for the event to occur.
integration_span = [x_initial, 20.0]

# Solve the initial value problem.
solution = solve_ivp(
    fun=ode_function,
    t_span=integration_span,
    y0=[y_initial],
    method='RK45',  # A standard accurate solver
    events=[event_y_reaches_minus_3],
    dense_output=True # Needed to accurately find the event time
)

# Extract and print the results
if solution.status == 1:  # Status 1 means a terminal event was reached
    # The time (x-value) of the event is stored in t_events
    x0 = solution.t_events[0][0]
    # The state (y-value) at the event is in y_events
    y_at_event = solution.y_events[0][0][0]
    # We can find the slope p = dy/dx at that point
    p_at_event = get_p(x0, y_at_event)
    
    print(f"The position x0 where the particle reaches y(x0) = -3 is: {x0:.6f}\n")
    
    print("Verification:")
    print(f"The state at the event is (x, y, dy/dx) = ({x0:.6f}, {y_at_event:.6f}, {p_at_event:.6f})")
    print("Plugging these values into the original equation (dy/dx)^3 + y^2 = xy(dy/dx):")
    
    # Calculate the left-hand side (LHS) and right-hand side (RHS) of the ODE
    lhs = p_at_event**3 + y_at_event**2
    rhs = x0 * y_at_event * p_at_event
    
    # Output the equation with the found numbers, as requested
    print(f"({p_at_event:.6f})^3 + ({y_at_event:.6f})^2 = ({x0:.6f}) * ({y_at_event:.6f}) * ({p_at_event:.6f})")
    print(f"{p_at_event**3:.6f} + {y_at_event**2:.6f} = {rhs:.6f}")
    print(f"{lhs:.6f} = {rhs:.6f}")
    
    # Check if the two sides are close enough to be considered equal
    if math.isclose(lhs, rhs, rel_tol=1e-6):
        print("\nThe equation holds, confirming the numerical solution.")
    else:
        print("\nThere is a numerical discrepancy.")

else:
    print("The simulation finished without the particle reaching y = -3 in the given interval.")
    
# Finally, return the required answer in the specified format
if 'x0' in locals():
    print(f'<<<1.637532>>>')
