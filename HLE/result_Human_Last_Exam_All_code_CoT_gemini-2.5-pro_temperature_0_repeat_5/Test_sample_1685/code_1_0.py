import numpy as np
from scipy.integrate import solve_ivp

def find_asymptotic_solution():
    """
    This function finds the approximate analytical solution for the given ODE in the large x regime.
    It first derives the form of the solution analytically and then solves for the unknown parameters
    by numerically integrating the full ODE.
    """

    # 1. Analytical derivation of constants A and q
    # From the dominant balance y''' ~ (y')^4, we find the form y(x) = C + A*(x0-x)^q
    # This gives q = 1/3 and A^3 = -30.
    q = 1/3
    A = (-30)**(1/3)

    # 2. Numerical simulation to find C and x0
    # Define the ODE system y_vec = [y, y', y'']
    def ode_system(x, y_vec):
        y, y_prime, y_pprime = y_vec
        
        # The tan(x) term can cause issues. We handle its singularities.
        # The main singularity is expected to be driven by the y'^4 term,
        # which likely occurs before the first singularity of tan(x) at pi/2.
        cos_x = np.cos(x)
        if abs(cos_x) < 1e-9: # Avoid division by zero for tan(x)
            tan_term = 0
        else:
            tan_val = np.tan(x)
            if abs(tan_val + 1) < 1e-9: # Avoid singularity in the ODE term
                # This is a singularity of the ODE itself. The solver should stop.
                # We return a large number to make the solver stop or slow down.
                return [y_prime, y_pprime, 1e10]
            tan_term = 1 / (tan_val + 1)

        y_ppprime = y**4 + y_prime**4 - y_pprime / (3*x**2 + 2) + tan_term
        return [y_prime, y_pprime, y_ppprime]

    # Initial conditions from the problem
    y0 = [0.00, 3.00, 2.00]
    
    # We expect a singularity, so we don't integrate too far.
    x_span = [0, 0.5]

    # We define an event to stop the integration when y' becomes large,
    # indicating we are approaching the singularity.
    def singularity_event(x, y_vec):
        return 1000 - y_vec[1] # Stop when y' reaches 1000
    singularity_event.terminal = True # Stop the integration when the event occurs
    singularity_event.direction = -1 # Event triggers when the value goes from positive to negative

    # Solve the ODE
    sol = solve_ivp(
        ode_system, 
        x_span, 
        y0, 
        dense_output=True, 
        events=singularity_event,
        rtol=1e-8, # Use high tolerance for accuracy
        atol=1e-8
    )

    # 3. Extract parameters and print the final expression
    if sol.status == 1 and sol.t_events[0].size > 0:
        # The event was triggered, we found the singularity
        x0_est = sol.t_events[0][0]
        # Get the solution vector [y, y', y''] at the singularity
        y_vec_at_x0 = sol.sol(x0_est)
        C_est = y_vec_at_x0[0]

        print("An approximate analytical expression for the solution y(x) in the large x regime (near the singularity) is:")
        # We format the output string with the calculated values rounded to two decimal places.
        # The final equation is y(x) = C + A * (x0 - x)^q
        # Note that A is negative, so we print it as C - |A| * ...
        print(f"y(x) = {C_est:.2f} + ({A:.2f}) * ({x0_est:.2f} - x)^({q:.2f})")
        print(f"y(x) = {C_est:.2f} - {abs(A):.2f} * ({x0_est:.2f} - x)^({q:.2f})")

    else:
        print("Numerical integration did not find a singularity in the specified range.")
        print("The asymptotic analysis might not be applicable or the singularity is further away.")

# Execute the function to find and print the solution
find_asymptotic_solution()