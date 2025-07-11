import numpy as np
from scipy.integrate import solve_ivp

# The differential equation is p^3 - x*y*p + y^2 = 0, where p = dy/dx.
# We need to define a function that returns p for a given (t, y) pair.
# Note: solve_ivp uses 't' as the independent variable, which is 'x' in our problem.
def ode_function(t, y):
    """
    Defines the differential equation dy/dt = p(t, y).
    """
    x = t
    # Coefficients of the cubic polynomial in p: p^3 + 0*p^2 - (x*y)*p + y^2 = 0
    coeffs = [1, 0, -x * y, y**2]
    
    # Find the roots of the cubic polynomial
    roots = np.roots(coeffs)
    
    # The roots can be complex. We are interested in the real root.
    # For the region of interest, there is only one real root.
    real_roots = roots[np.isreal(roots)].real
    
    # We expect only one real root based on the discriminant analysis.
    # If there are multiple, this problem would require careful selection
    # to maintain continuity, but here we can just take the first (and only).
    p = real_roots[0]
    
    return p

def event_y_equals_minus_3(t, y):
    """
    Event function to find when y = -3.
    The solver will find the root of this function.
    """
    return y[0] + 3

# We want the solver to stop when it finds the root, so set terminal=True
event_y_equals_minus_3.terminal = True
# The event function is zero when y is going from >-3 to <-3, which is a negative direction.
event_y_equals_minus_3.direction = -1

# Initial conditions
x0 = 0
y0 = -1

# We need to integrate over a sufficiently large span for x to find the solution
t_span = (x0, 10.0) # Integrate from x=0 up to x=10
y_init = [y0]

# Solve the ODE
# We use a high accuracy setting to get a precise result.
sol = solve_ivp(
    fun=ode_function,
    t_span=t_span,
    y0=y_init,
    events=event_y_equals_minus_3,
    dense_output=True,
    rtol=1e-8,
    atol=1e-8
)

# The time at which the event occurred is stored in sol.t_events
if sol.status == 1: # Status 1 means a terminal event was triggered
    x_solution = sol.t_events[0][0]
    y_solution = sol.y_events[0][0]
    print(f"The particle reaches the vertical coordinate y = {y_solution:.1f} at position x = {x_solution:.1f}.")
    print("The final equation for the particle's position is: ")
    print(f"x0 = {x_solution:.1f}")

else:
    print("The solution was not found in the integration span.")

# Let's verify the solution by finding the slope p at the solution point (x,y) = (4.0, -3.0)
x_final = 4.0
y_final = -3.0
# p^3 + y^2 = xy*p => p^3 - (x*y)*p + y^2 = 0
# p^3 - (4.0 * -3.0)*p + (-3.0)^2 = 0
# p^3 + 12*p + 9 = 0
p_coeffs = [1, 0, 12, 9]
p_roots = np.roots(p_coeffs)
p_final_real = p_roots[np.isreal(p_roots)].real[0]
print(f"At the point (x, y) = ({x_final:.1f}, {y_final:.1f}), the equation for the slope 'p' is:")
print(f"p^3 + {int(x_final * -y_final)}*p + {int(y_final**2)} = 0")
print(f"The real solution for the slope at this point is p ≈ {p_final_real:.4f}")
