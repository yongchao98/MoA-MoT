import numpy as np
from scipy.integrate import solve_ivp

def solve_ode_system():
    """
    Solves the given system of ODEs to find the time t when b(t) = 0.5.
    """
    # Define the system of differential equations
    def system(t, y, A, k):
        """
        Defines the system of ODEs:
        a'(t) = -0.5*a^2 - A*b^2 + k*(b-1)
        b'(t) = -a*b
        y is a list or array [a(t), b(t)]
        """
        a, b = y
        da_dt = -0.5 * a**2 - A * b**2 + k * (b - 1)
        db_dt = -a * b
        return [da_dt, db_dt]

    # Parameters
    A = 1
    k = 5

    # Initial conditions
    a0 = 0.1
    b0 = 2.0
    y0 = [a0, b0]

    # Time span for integration.
    # Choose a sufficiently long time to ensure the event is found.
    t_span = [0, 20]

    # Define the event function to find when b(t) = 0.5
    def event_b_equals_0_5(t, y):
        """
        Event function that returns 0 when b(t) = 0.5.
        The solver will find the root of this function.
        """
        return y[1] - 0.5

    # This attribute tells the solver to stop integration when the event occurs
    event_b_equals_0_5.terminal = True
    # This attribute helps the solver find the event more efficiently.
    # Since b(t) is expected to decrease, the event function value will
    # be decreasing (going from positive to negative).
    event_b_equals_0_5.direction = -1

    # Solve the ODE using scipy's solve_ivp
    sol = solve_ivp(
        fun=lambda t, y: system(t, y, A, k),
        t_span=t_span,
        y0=y0,
        events=event_b_equals_0_5,
        dense_output=True  # Allows for accurate interpolation
    )

    # Output the results
    print("Solving the system of differential equations:")
    print(f"a'(t) = -0.5*a(t)^2 - {A}*b(t)^2 + {k}*(b(t)-1)")
    print(f"b'(t) = -a(t)*b(t)")
    print(f"\nWith initial conditions a(0) = {a0} and b(0) = {b0}.")
    
    # Check if the event was found and print the time
    if sol.status == 1 and sol.t_events[0].size > 0:
        t_event = sol.t_events[0][0]
        print(f"\nThe time t at which b(t) = 0.5 is estimated to be: {t_event:.4f}")
    else:
        print(f"\nThe condition b(t) = 0.5 was not met in the time interval {t_span}.")
        print("This could mean no such t exists or a larger time interval is needed.")

if __name__ == '__main__':
    solve_ode_system()
<<<A>>>