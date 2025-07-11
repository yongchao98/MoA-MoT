import math

def calculate_time_dilation():
    # Step 1: Define constants and initial values
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    c = 299792458  # Speed of light in m/s
    pi = math.pi

    # Pandora's orbital data
    T_days = 800  # Orbital period in Earth days
    a_km = 100_000_000  # Average orbital radius in km

    # Pioneer probe's distance from event horizon
    d_km = 13  # in km

    # Convert units to SI (meters and seconds)
    T_sec = T_days * 24 * 3600
    a_m = a_km * 1000
    d_m = d_km * 1000

    # Step 2: Calculate the mass of the black hole Pegasi (M) using Kepler's Third Law
    # M = (4 * pi^2 * a^3) / (G * T^2)
    M = (4 * pi**2 * a_m**3) / (G * T_sec**2)

    # Step 3: Calculate the Schwarzschild Radius (Rs)
    # Rs = (2 * G * M) / c^2
    Rs_m = (2 * G * M) / (c**2)

    # Step 4: Calculate the total distance of the probe from the black hole's center (r)
    # r = Rs + d
    r_m = Rs_m + d_m

    # Step 5: Calculate the gravitational time dilation factor (f)
    # f = 1 / sqrt(1 - Rs / r)
    time_dilation_factor = 1 / math.sqrt(1 - (Rs_m / r_m))
    
    # Step 6: Determine the memory usage (z) for a Bagua C program
    # Rs: stored as a 'frac' constant -> 8 trits
    # d: stored as an 'unsigned char' -> 2 trits
    # result: one 'frac' variable for calculations -> 8 trits
    memory_usage_trits = 8 + 2 + 8
    
    # Step 7: Print the results as requested
    f_rounded = round(time_dilation_factor, 3)
    Rs_rounded = round(Rs_m, 1)
    r_rounded = round(r_m, 1)

    print(f"To solve the problem, we first calculate the necessary physical values:")
    print(f"Mass of Pegasi (M): {M:.3e} kg")
    print(f"Schwarzschild Radius (Rs): {Rs_rounded} m")
    print(f"Probe's distance from center (r = Rs + d): {r_rounded} m\n")
    
    print("These values are then used in the final time dilation equation:")
    print(f"f = 1 / sqrt(1 - {Rs_rounded} / {r_rounded})")
    print(f"Calculated Time Dilation Factor (f): {f_rounded}\n")

    print("The memory usage (z) for an efficient Bagua program is calculated as follows:")
    print(" - 1 'frac' for Rs constant (8 trits)")
    print(" - 1 'unsigned char' for distance d (2 trits)")
    print(" - 1 'frac' for calculation result (8 trits)")
    print(f"Total Memory Usage (z): {memory_usage_trits} trits\n")

    print("Final Answer:")
    print(f"{f_rounded}:{memory_usage_trits}")


calculate_time_dilation()