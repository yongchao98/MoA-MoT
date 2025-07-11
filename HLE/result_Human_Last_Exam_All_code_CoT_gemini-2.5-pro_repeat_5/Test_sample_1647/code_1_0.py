import math

def calculate_time_dilation_and_memory():
    """
    Calculates the time dilation factor and memory usage for the Pioneer probe problem.
    """
    # Step 1: Define constants and convert to SI units
    # Pandora's orbital period T = 800 days
    T_seconds = 800 * 24 * 3600

    # Pandora's average orbital radius a = 100,000,000 km
    a_meters = 100e6 * 1000

    # Pioneer's distance from event horizon d = 13 km
    d_meters = 13 * 1000

    # Physical constants
    c_mps = 299792458  # Speed of light in m/s

    # Step 2: Calculate the Schwarzschild Radius (Rs) of the black hole Pegasi
    # Using the formula derived from Kepler's law: Rs = (8 * pi^2 * a^3) / (c^2 * T^2)
    numerator_rs = 8 * (math.pi**2) * (a_meters**3)
    denominator_rs = (c_mps**2) * (T_seconds**2)
    Rs_meters = numerator_rs / denominator_rs

    # Step 3: Calculate the time dilation factor f
    # The most memory-efficient C program would use the Taylor approximation f ≈ 1 + (Rs/d)/2.
    # We will calculate the precise value here for accuracy.
    # f = sqrt(1 + Rs/d)
    time_dilation_factor = math.sqrt(1 + Rs_meters / d_meters)

    # Round f to 3 decimal places as required
    f_rounded = round(time_dilation_factor, 3)

    # Step 4: Calculate the memory usage (z) for the most efficient Bagua C program
    # The program would be:
    #   frac term = ...; // Pre-calculated (Rs/d)/2
    #   frac f = 1 + term;
    # This requires 2 variables of type 'frac'.
    
    # Size of Bagua data types in trits
    trit_size_char = 2
    trit_size_wchar = 4

    # Size of one 'frac' struct: signed char (n) + unsigned wchar (d) + signed char (e)
    trit_size_frac = trit_size_char + trit_size_wchar + trit_size_char
    
    # Total memory usage for 2 'frac' variables
    num_variables = 2
    z_memory_usage = num_variables * trit_size_frac
    
    # Step 5: Print the final result in the format f:z
    # The final equation's numbers are the calculated factor and memory usage.
    print(f"{f_rounded}:{z_memory_usage}")

calculate_time_dilation_and_memory()