import math

def calculate_critical_speed():
    """
    This script calculates the critical speed for a vehicle using the linear
    single-track (bicycle) model.

    Plan:
    1. Define the parameters of the vehicle model.
    2. Check the condition for oversteering (a*c_f > b*c_r).
    3. If the vehicle is oversteering, calculate the critical speed using the
       derived formula:
       v_crit = sqrt(((a+b)^2 * c_f * c_r) / (m * (a*c_f - b*c_r)))
    4. Print the formula with the parameter values substituted in.
    5. Print the final calculated critical speed in both m/s and km/h.
    6. If the vehicle is not oversteering (i.e., understeering or neutral steer),
       print a message indicating that it is stable according to this model.
    """
    
    # --- Vehicle Parameters ---
    # You can change these values to match your specific vehicle
    a = 1.2      # distance from CG to front axle (m)
    b = 1.6      # distance from CG to rear axle (m)
    c_f = 100000 # cornering stiffness of the front axle (N/rad)
    c_r = 50000  # cornering stiffness of the rear axle (N/rad)
    m = 1500     # vehicle mass (kg)
    # Note: Moment of inertia (I) is not needed for the critical speed formula.
    
    print("--- Vehicle Parameters ---")
    print(f"a (front CG distance) = {a} m")
    print(f"b (rear CG distance)  = {b} m")
    print(f"c_f (front stiffness) = {c_f} N/rad")
    print(f"c_r (rear stiffness)  = {c_r} N/rad")
    print(f"m (mass)              = {m} kg")
    print("-" * 30)

    # Step 1: Check for oversteering condition
    oversteer_term = a * c_f - b * c_r
    
    print("Checking steering characteristic (a*c_f - b*c_r):")
    print(f"({a} * {c_f}) - ({b} * {c_r}) = {a*c_f} - {b*c_r} = {oversteer_term}")

    # Step 2: Calculate critical speed based on the steering characteristic
    if oversteer_term > 0:
        print("\nThe vehicle is OVERSTEERING. A critical speed exists.")
        print("Calculating the speed at which it becomes unstable...")
        
        # Calculate the components of the formula
        wheelbase = a + b
        numerator = (wheelbase**2) * c_f * c_r
        denominator = m * oversteer_term
        
        # Calculate the critical speed
        v_crit_sq = numerator / denominator
        v_crit_ms = math.sqrt(v_crit_sq) # speed in m/s
        v_crit_kmh = v_crit_ms * 3.6    # speed in km/h
        
        # Print the final equation with values substituted
        print("\n--- Calculation using the formula ---")
        print(f"v_crit = sqrt( ( (a+b)^2 * c_f * c_r ) / ( m * (a*c_f - b*c_r) ) )")
        print(f"v_crit = sqrt( ( ({a} + {b})^2 * {c_f} * {c_r} ) / ( {m} * ({oversteer_term}) ) )")
        print(f"v_crit = sqrt( ( {numerator:.0f} ) / ( {denominator:.0f} ) )")
        print(f"v_crit = sqrt( {v_crit_sq:.2f} )")
        
        # Print the final result
        print("\n--- Final Critical Speed ---")
        print(f"v_crit = {v_crit_ms:.2f} m/s")
        print(f"v_crit = {v_crit_kmh:.2f} km/h")

    else:
        # Handle understeer or neutral steer cases
        if oversteer_term < 0:
            print("\nThe vehicle is UNDERSTEERING. According to this linear model, it is stable at all speeds.")
        else: # oversteer_term == 0
            print("\nThe vehicle is NEUTRAL STEER. According to this linear model, it is stable at all speeds.")

# Execute the calculation
calculate_critical_speed()