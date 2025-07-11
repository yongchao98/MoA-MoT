import numpy as np
from scipy.integrate import solve_ivp

def main():
    """
    This script solves for the instantaneous rate of change of a pendulum's angle.
    """
    # Define the differential equation dy/dt = f(t, y)
    def pendulum_eq(t, y):
        """
        Defines the differential equation for the pendulum's motion.
        dy/dt = sin(y(t) - t) - 4
        """
        # The input y is an array, so we access its first element.
        return np.sin(y[0] - t) - 4

    # Set the initial condition: y(0) = 0
    y0 = [0]

    # Set the time at which we want to find the rate of change
    t_final = np.pi / 6

    # Set the time span for the numerical integration
    t_span = [0, t_final]

    # Use the numerical solver to find the value of y at t_final.
    # We request the solution specifically at the endpoint t_final.
    solution = solve_ivp(pendulum_eq, t_span, y0, t_eval=[t_final])

    # Extract the calculated value of y(π/6) from the solution object.
    # solution.y is a 2D array of shape (number_of_equations, number_of_time_points)
    y_final = solution.y[0, 0]

    # Now, calculate the instantaneous rate of change dy/dt at t = π/6
    # using the given formula: dy/dt = sin(y(t) - t) - 4
    rate_of_change = np.sin(y_final - t_final) - 4

    # Print the steps of the calculation as requested
    print("The equation for the rate of change is: dy/dt = sin(y(t) - t) - 4")
    print(f"We need to find the value of dy/dt at t = π/6.")
    print(f"\nStep 1: Find the value of y(t) at t = π/6.")
    print(f"Using a numerical solver with the initial condition y(0) = 0:")
    print(f"The value of t is π/6 ≈ {t_final:.8f}")
    print(f"The calculated value of y at this time is y(π/6) ≈ {y_final:.8f}")

    print(f"\nStep 2: Substitute these values into the equation for dy/dt.")
    # Show each number in the final equation
    sin_argument = y_final - t_final
    sin_value = np.sin(sin_argument)
    
    print(f"dy/dt = sin(y(π/6) - π/6) - 4")
    print(f"dy/dt = sin({y_final:.8f} - {t_final:.8f}) - 4")
    print(f"dy/dt = sin({sin_argument:.8f}) - 4")
    print(f"dy/dt = {sin_value:.8f} - 4")
    print(f"\nThe final calculated instantaneous rate of change is: {rate_of_change:.8f}")

if __name__ == "__main__":
    main()