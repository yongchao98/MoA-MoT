import math

def solve_time_dilation():
    """
    Calculates the gravitational time dilation factor 'f' and the memory usage 'z'
    for the most efficient Bagua C program.
    """

    # 1. Define physical constants and problem parameters in SI units.
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    c = 299792458    # Speed of light in m/s
    pi = math.pi
    
    # Pandora's orbital data
    pandora_T_days = 800
    pandora_a_km = 100_000_000
    
    # Probe's distance from event horizon
    probe_d_km = 13

    # Convert parameters to SI units (meters, kilograms, seconds)
    T_sec = pandora_T_days * 24 * 3600
    a_m = pandora_a_km * 1000
    d_m = probe_d_km * 1000

    # 2. Calculate the mass of the black hole Pegasi (M) using Kepler's Third Law.
    M_kg = (4 * pi**2 * a_m**3) / (G * T_sec**2)

    # 3. Calculate the Schwarzschild radius (Rs) of the black hole.
    Rs_m = (2 * G * M_kg) / (c**2)

    # 4. Calculate the time dilation factor 'f'.
    # 'r' is the probe's distance from the center of the black hole.
    r_m = Rs_m + d_m
    
    # Use the exact formula for gravitational time dilation.
    ratio = Rs_m / r_m
    if ratio >= 1:
        # This case corresponds to being inside the event horizon.
        f = float('inf')
    else:
        f = 1 / math.sqrt(1 - ratio)

    # Round the calculated factor 'f' to three decimal places.
    f_rounded = round(f, 3)

    # 5. Determine the memory usage 'z'.
    # The most memory-efficient C program on Bagua would use a single 'frac' variable
    # to perform a Taylor series approximation of the sqrt function.
    # The size of a `frac` variable is specified as 8 trits.
    z_trits = 8

    # 6. Print the final answer. The 'final equation' is the output format f:z.
    # We output each number, f and z, required for this result.
    print(f"{f_rounded}:{z_trits}")

solve_time_dilation()