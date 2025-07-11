import numpy as np
from scipy.integrate import solve_ivp

def solve_and_estimate_area():
    """
    Solves the system of ODEs for a grid of initial conditions and estimates
    the area of the set Omega where solutions blow up as specified.
    """
    
    # 1. Define the system of ODEs
    def ode_system(t, y):
        a, b = y
        a_prime = -b * a
        b_prime = -b**2 / 2 - a**2 + 6 * (a - 1)
        return [a_prime, b_prime]

    # 2. Set up the grid of initial conditions
    a_range = [-1, 1]
    b_range = [2, 3]
    # Using a high resolution grid for better accuracy
    n_a = 201 
    n_b = 101 
    a_vals = np.linspace(a_range[0], a_range[1], n_a)
    b_vals = np.linspace(b_range[0], b_range[1], n_b)
    
    total_points = n_a * n_b
    blow_up_count = 0

    # 3. Set up the numerical integration parameters
    t_max = 20.0
    t_span = [0, t_max]
    blow_up_thresh = 1000.0

    # Define an event to stop integration if solution gets too large
    def blow_up_event(t, y):
        # Event triggers when |a| or |b| exceeds the threshold
        return np.max(np.abs(y)) - blow_up_thresh
    blow_up_event.terminal = True  # Stop integration when event occurs
    
    # 4. Iterate over initial conditions and solve the ODE
    for a0 in a_vals:
        for b0 in b_vals:
            y0 = [a0, b0]
            
            sol = solve_ivp(
                ode_system, 
                t_span, 
                y0, 
                events=blow_up_event,
                dense_output=True # dense_output is good practice with events
            )
            
            # 5. Check if the blow-up condition is met
            # The solver stops successfully when an event is triggered (sol.status == 1)
            # We check if the event array is not empty.
            if sol.status == 1 and len(sol.t_events[0]) > 0:
                # Get the state at the time of the event
                a_final, b_final = sol.y_events[0][0]
                
                # The desired blow-up is a -> +inf and b -> -inf.
                # This is identified by a > 0 and b < 0 at the time of blow-up.
                if a_final > 0 and b_final < 0:
                    blow_up_count += 1
                    
    # 6. Estimate the measure of the set Omega
    # Total area of the initial parameter rectangle
    total_area = (a_range[1] - a_range[0]) * (b_range[1] - b_range[0])
    
    # Estimated area is the fraction of blow-up points times the total area
    m_omega_estimated = (blow_up_count / total_points) * total_area
    
    print(f"Initial domain R is [a,b] in [-1, 1] x [2, 3], with total area = {total_area}")
    print(f"Grid size used: {n_a} (for a) x {n_b} (for b) = {total_points} total initial points.")
    print(f"Number of initial points in R that lead to the desired blow-up: {blow_up_count}")
    print(f"The estimated measure (area) of the set Omega is: {m_omega_estimated}")

# Execute the function to get the answer
solve_and_estimate_area()