import math

# Define the parameters of the linear single-track model
# a: distance from CG to front axle (m)
# b: distance from CG to rear axle (m)
# cf: cornering stiffness of the front axle (N/rad)
# cr: cornering stiffness of the rear axle (N/rad)
# m: vehicle mass (kg)

# Example values for an oversteering vehicle
a = 1.2   # m
b = 1.5   # m
cf = 80000 # N/rad
cr = 60000 # N/rad
m = 1500  # kg

# --- Step 1: Check for oversteering condition ---
# The critical speed exists only for oversteering vehicles, where a*cf > b*cr.
is_oversteering = (a * cf) > (b * cr)

if not is_oversteering:
    print("The vehicle is not oversteering based on the given parameters.")
    print(f"Condition check: a*cf = {a*cf:.2f}, b*cr = {b*cr:.2f}")
    print("An understeering or neutral steer vehicle does not have a critical speed.")
else:
    # --- Step 2: Calculate the critical speed ---
    print("The vehicle is oversteering. Calculating the critical speed.")
    print(f"Condition check: a*cf = {a * cf:.0f} > b*cr = {b * cr:.0f}\n")
    
    # Numerator of the term inside the square root
    numerator = ((a + b)**2) * cf * cr
    
    # Denominator of the term inside the square root
    denominator = m * (a * cf - b * cr)
    
    # Calculate the square of the critical speed
    v_crit_squared = numerator / denominator
    
    # Calculate the critical speed
    v_crit = math.sqrt(v_crit_squared)
    
    # --- Step 3: Display the derivation and result ---
    print("The formula for critical speed (v_crit) is:")
    print("v_crit = sqrt( ((a + b)^2 * cf * cr) / (m * (a * cf - b * cr)) )\n")
    
    print("Substituting the given values into the equation:")
    print(f"v_crit = sqrt( (({a} + {b})^2 * {cf} * {cr}) / ({m} * ({a} * {cf} - {b} * {cr})) )")
    print(f"v_crit = sqrt( ({numerator:.2e}) / ({denominator:.2e}) )")
    print(f"v_crit = sqrt( {v_crit_squared:.2f} )")
    print(f"v_crit = {v_crit:.2f} m/s")
    
    # Also print in km/h for better intuition
    v_crit_kmh = v_crit * 3.6
    print(f"v_crit = {v_crit_kmh:.2f} km/h")
    
    # Final answer for the user
    # print(f"\n<<< {v_crit:.2f} >>>") # This is a placeholder for the final format
    
# The problem asks to output the final result in a specific format,
# which is the numerical value of the critical speed in m/s.
final_answer = v_crit if is_oversteering else "N/A"
# The final answer is wrapped in <<<>>>
# print(f'<<<{v_crit if is_oversteering else "N/A"}>>>')
# For this task, I will output the final result as a number
if is_oversteering:
    print(f"\nFinal Answer (in m/s):")
    print(f'<<<{v_crit:.2f}>>>')
