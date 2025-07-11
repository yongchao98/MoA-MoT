import numpy as np
from scipy.integrate import solve_ivp
import warnings

# We can ignore warnings about overflow as they are expected when the solution blows up.
warnings.filterwarnings("ignore")

def find_approximation():
    """
    This script finds the approximate analytical expression for the given ODE
    by numerically calculating the point of singularity, x0.
    """
    
    # Define the system of first-order ODEs
    # z = [y, y', y'']
    def ode_system(x, z):
        y, yp, ypp = z
        
        # This term has poles, but the solution is expected to blow up before the first pole.
        tan_term = 1.0 / (np.tan(x) + 1.0)
        
        yppp = y**4 + yp**4 - ypp / (3*x**2 + 2) + tan_term
        return [yp, ypp, yppp]

    # Initial conditions
    z0 = [0.00, 3.00, 2.00]

    # Time span for integration. The first pole of the tan term is at 3*pi/4 (~2.356).
    # We expect the blow-up to happen before that.
    t_span = [0, 2.3]

    # Event function to detect the blow-up singularity.
    # We define blow-up as y reaching a very large value.
    def blowup_event(x, z):
        return z[0] - 1e6
    blowup_event.terminal = True # Stop integration when event occurs
    blowup_event.direction = 1   # Event triggers when z[0] is increasing

    # Solve the ODE
    sol = solve_ivp(ode_system, t_span, z0, events=blowup_event, dense_output=True)

    # Extract the results
    if sol.status == 1 and sol.t_events[0].size > 0:
        x0 = sol.t_events[0][0]
        
        # From the asymptotic analysis, the constant C is 6^(1/3)
        C = 6**(1/3)

        # Round the constants to two decimal places
        C_rounded = round(C, 2)
        x0_rounded = round(x0, 2)

        print("Based on the analysis, the approximate expression is of the form: y(x) = C / (x - x0)")
        print("The values of the constants are:")
        print(f"C = {C_rounded}")
        print(f"x0 = {x0_rounded}")
    else:
        print("The solution did not blow up in the integration interval.")
        print("An approximate expression of the assumed form could not be determined.")

# Run the function to find and print the result
find_approximation()