import numpy as np
from scipy.integrate import solve_ivp

def solve_trajectory():
    """
    Solves the given ODE numerically to find the x-coordinate where y = -3.
    """

    # Define the function f(x, y) which returns dy/dx.
    # To do this, we must solve the cubic equation for p = dy/dx at each step.
    # The equation is: p^3 - (x*y)*p + y^2 = 0
    def ode_func(x, y_arr):
        y = y_arr[0]
        # Coefficients of the cubic polynomial in p (dy/dx)
        coeffs = [1, 0, -x * y, y**2]
        
        # Find the roots of the cubic equation
        roots = np.roots(coeffs)
        
        # There will always be at least one real root. We select it.
        # A robust way is to find the root with the smallest imaginary part.
        real_root = roots[np.isreal(roots)].real
        if len(real_root) == 0: # Should not happen for a cubic, but for safety
            real_root = roots[np.argmin(np.abs(roots.imag))].real
        else:
            # If multiple real roots, we need to choose the one that continues the trajectory.
            # At (0, -1), p must be -1. We assume the solution is continuous,
            # so we'll pick the root closest to the previous step's slope if needed.
            # For this problem, it turns out there's only one real root along the path.
            real_root = real_root[0]
            
        return real_root

    # Define an event to find the point where y(x) = -3.
    # The event function should return 0 at the event.
    def reach_y_minus_3(x, y):
        return y[0] + 3

    # At t=0, y=-1, so y+3=2. The slope is initially negative, so y decreases.
    # The event function will go from positive to negative.
    reach_y_minus_3.direction = -1
    # We want to stop the integration when this happens.
    reach_y_minus_3.terminal = True

    # Initial conditions
    y0 = [-1.0]
    x_span = (0, 5)  # Integrate over a reasonable interval for x

    # Solve the ODE
    sol = solve_ivp(
        fun=ode_func,
        t_span=x_span,
        y0=y0,
        method='RK45',
        events=reach_y_minus_3,
        dense_output=True,
        rtol=1e-8,
        atol=1e-8
    )

    # Check if the event was triggered
    if sol.t_events[0].size > 0:
        x0 = sol.t_events[0][0]
        y0_final = -3.0
        
        # Calculate the slope p0 = dy/dx at the final point
        p0 = ode_func(x0, [y0_final])
        
        print(f"The particle reaches the vertical coordinate y = {y0_final} at position x0 = {x0:.8f}")
        print("\nThe governing equation (dy/dx)^3 + y^2 = xy(dy/dx) is satisfied at this point:")
        
        # Print the equation with the found values
        print(f"({p0:.8f})^3 + ({y0_final:.1f})^2 = ({x0:.8f}) * ({y0_final:.1f}) * ({p0:.8f})")
        
        lhs = p0**3 + y0_final**2
        rhs = x0 * y0_final * p0
        print(f"{lhs:.8f} = {rhs:.8f}")
        
        print(f"<<<{x0}>>>")
    else:
        print(f"The trajectory did not reach y = -3 within the x interval {x_span}.")

solve_trajectory()