import numpy as np
from scipy.integrate import solve_ivp

def get_solution_behavior():
    """
    Numerically estimates the measure of the set Omega by sampling initial conditions.
    """
    
    # Define the system of ODEs y = [a, b]
    def ode_system(t, y):
        a, b = y
        try:
            # The exponential term can grow very quickly.
            et = np.exp(t)
        except OverflowError:
            et = np.inf
        
        da_dt = -b * a
        db_dt = -0.5 * b**2 - et * a**2 - a
        return [da_dt, db_dt]

    # Event function to stop integration if variables get too large (potential blow-up)
    def blow_up_event(t, y):
        # Stop if |a| or |b| exceeds a large threshold
        return max(abs(y[0]), abs(y[1])) - 1e6
    blow_up_event.terminal = True
    blow_up_event.direction = 1

    # As per our analysis, we only need to test a0 > 0.
    # The domain for these initial conditions is (0, 1] x [10, 20].
    candidate_area = (1.0 - 0.0) * (20.0 - 10.0)

    # Create a grid of initial conditions in the candidate domain.
    # We start a0 from a small positive number to avoid a0=0.
    a0_vals = np.linspace(0.1, 1.0, 10)
    b0_vals = np.linspace(10.0, 20.0, 11)
    
    total_points_tested = len(a0_vals) * len(b0_vals)
    blow_up_count = 0
    t_max = 50.0  # Maximum integration time, blow-up should happen before this.

    # Iterate over the grid of initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            y0 = [a0, b0]
            
            # Solve the ODE system
            sol = solve_ivp(
                ode_system, 
                [0, t_max], 
                y0, 
                events=blow_up_event,
                method='Radau'  # A good choice for potentially stiff problems
            )

            # Check if a blow-up event occurred
            if sol.status == 1 and sol.t_events[0].size > 0:
                # Event was triggered. Check the state at the event time.
                a_final, b_final = sol.y_events[0][-1]
                
                # Our specific blow-up condition: a -> +inf, b -> -inf
                if a_final > 1e5 and b_final < -1e5:
                    blow_up_count += 1
    
    # Estimate the measure of Omega
    if total_points_tested > 0:
        blow_up_fraction = blow_up_count / total_points_tested
    else:
        blow_up_fraction = 0
        
    estimated_area = blow_up_fraction * candidate_area

    # Print the results of the estimation calculation
    print("Numerical Estimation of m(Omega):")
    print(f"Candidate region for initial conditions (a0, b0): (0, 1] x [10, 20]")
    print(f"Area of candidate region: {candidate_area}")
    print(f"Number of initial points tested: {total_points_tested}")
    print(f"Number of points leading to blow-up: {blow_up_count}")
    print("\nFinal calculation for the estimated measure:")
    # The final equation as requested.
    print(f"m(Omega) = (Number of blow-ups / Total points tested) * Area of candidate region")
    print(f"m(Omega) = ({blow_up_count} / {total_points_tested}) * {candidate_area}")
    print(f"m(Omega) = {blow_up_fraction:.4f} * {candidate_area}")
    print(f"Estimated m(Omega) = {estimated_area:.4f}")

# Execute the function to get the answer
get_solution_behavior()