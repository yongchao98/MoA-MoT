import math

# --- Vehicle Parameters ---
# a: distance from CG to front axle (m)
# b: distance from CG to rear axle (m)
# cf: cornering stiffness of the front axle (N/rad)
# cr: cornering stiffness of the rear axle (N/rad)
# m: vehicle mass (kg)
# I: vehicle moment of inertia (kg*m^2) - Note: not needed for the final formula

# Example values for an oversteering vehicle
a = 1.5   # meters
b = 1.3   # meters
cf = 90000  # N/rad
cr = 70000  # N/rad
m = 1500  # kg
I = 2800  # kg*m^2

print("Deriving the critical speed for an oversteering vehicle.\n")
print("--- Given Parameters ---")
print(f"a  = {a} m")
print(f"b  = {b} m")
print(f"cf = {cf} N/rad")
print(f"cr = {cr} N/rad")
print(f"m  = {m} kg")
print("-" * 25 + "\n")

# --- Step 1: Check for Oversteering Condition ---
# For a critical speed to exist, the vehicle must be oversteering.
# The condition for oversteering is: a*cf > b*cr
oversteer_check_value = a * cf - b * cr

print("--- Step 1: Check for Oversteering Condition ---")
print("Condition for oversteer: a*cf - b*cr > 0")
print(f"Calculating a*cf - b*cr = {a} * {cf} - {b} * {cr} = {a*cf} - {b*cr} = {oversteer_check_value}")

if oversteer_check_value <= 0:
    print("\nThe vehicle is not oversteering. It is understeering or neutral steering.")
    print("A real critical speed does not exist (theoretically, it is infinite).")
else:
    print("The condition is met. The vehicle is oversteering.\n")

    # --- Step 2: Calculate Critical Speed using the Formula ---
    print("--- Step 2: Calculate Critical Speed ---")
    print("The formula for critical speed (v_crit) is:")
    print("v_crit = sqrt( ( (a+b)^2 * cf * cr ) / ( m * (a*cf - b*cr) ) )\n")
    
    # Numerator calculation
    L = a + b
    numerator = (L**2) * cf * cr
    
    # Denominator calculation
    denominator = m * oversteer_check_value

    print("Substituting the values into the formula:")
    print(f"v_crit = sqrt( ( ({a} + {b})^2 * {cf} * {cr} ) / ( {m} * ({a} * {cf} - {b} * {cr}) ) )")
    print(f"v_crit = sqrt( ( {L}^2 * {cf * cr} ) / ( {m} * ({a*cf} - {b*cr}) ) )")
    print(f"v_crit = sqrt( ( {L**2:.2f} * {cf * cr} ) / ( {m} * {oversteer_check_value} ) )")
    print(f"v_crit = sqrt( {numerator} / {denominator} )")

    # Final calculation
    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)

    print(f"v_crit = sqrt( {v_crit_sq:.4f} )")
    print(f"v_crit = {v_crit:.4f} m/s")

    # Convert to km/h for context
    v_crit_kmh = v_crit * 3.6
    print(f"\nIn other units: {v_crit_kmh:.2f} km/h")
    
    # Final answer in m/s as requested
    print(f"\n<<<v_crit = {v_crit}>>>")
