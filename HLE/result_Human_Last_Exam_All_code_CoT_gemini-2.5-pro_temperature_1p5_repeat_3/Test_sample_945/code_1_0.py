import math

def calculate_critical_speed():
    """
    Calculates and prints the critical speed for an oversteering vehicle based on
    the linear single-track model parameters.
    """
    # --- Parameters for a sample oversteering vehicle ---
    # a: distance from Center of Gravity (CG) to front axle (m)
    a = 1.2
    # b: distance from Center of Gravity (CG) to rear axle (m)
    b = 1.6
    # m: vehicle mass (kg)
    m = 1500
    # cf: front axle cornering stiffness (N/rad)
    cf = 80000
    # cr: rear axle cornering stiffness (N/rad)
    cr = 50000

    # The condition for an oversteering vehicle is: a * cf > b * cr
    oversteer_term = a * cf - b * cr

    print("--- Vehicle Parameter Analysis ---")
    print(f"Vehicle parameters:\na = {a} m\nb = {b} m\nm = {m} kg\ncf = {cf} N/rad\ncr = {cr} N/rad\n")

    print("--- Oversteering Condition Check ---")
    print(f"The oversteer condition is (a * cf - b * cr) > 0.")
    print(f"Value = ({a} * {cf}) - ({b} * {cr}) = {a * cf} - {b * cr} = {oversteer_term:.2f}")

    if oversteer_term <= 0:
        print("\nResult: The provided parameters do not describe an oversteering vehicle.")
        print("A critical speed does not exist under these conditions.")
        return

    print("Result: The condition is met, so the vehicle is oversteering.\n")

    # --- Critical Speed Calculation ---
    print("--- Critical Speed Derivation ---")
    print("The formula for critical speed (vcrit) is:")
    print("vcrit = sqrt( (cf * cr * (a + b)^2) / (m * (a * cf - b * cr)) )\n")

    # Calculate numerator and denominator for clarity
    numerator = cf * cr * (a + b)**2
    denominator = m * oversteer_term
    
    print("--- Calculation with Given Parameters ---")
    print("Plugging in the numbers into the formula:")
    print(f"vcrit = sqrt( ({cf} * {cr} * ({a} + {b})^2) / ({m} * ({a} * {cf} - {b} * {cr})) )")
    print(f"vcrit = sqrt( ({cf} * {cr} * {round((a + b)**2, 2)}) / ({m} * {oversteer_term}) )")
    print(f"vcrit = sqrt( {numerator:.2e} / {denominator:.2e} )")
    
    # Calculate final result
    vcrit_squared = numerator / denominator
    vcrit_ms = math.sqrt(vcrit_squared)
    vcrit_kmh = vcrit_ms * 3.6

    print(f"vcrit = sqrt({vcrit_squared:.2f}) m^2/s^2\n")

    print("--- Final Result ---")
    print(f"The critical speed is: {vcrit_ms:.2f} m/s")
    print(f"Which is equivalent to: {vcrit_kmh:.2f} km/h")


if __name__ == '__main__':
    calculate_critical_speed()