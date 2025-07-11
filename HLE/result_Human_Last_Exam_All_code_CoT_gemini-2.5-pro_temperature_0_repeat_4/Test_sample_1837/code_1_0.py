import math

def solve_star_measurement():
    """
    This script checks the consistency of stellar measurements using Planck's Law.
    It determines if the measurements are valid, and if not, identifies the most
    likely incorrect quantity and calculates its corrected value.
    """
    # --- Physical Constants ---
    h = 6.62607015e-34  # Planck constant (J*s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23   # Boltzmann constant (J/K)
    b = 2.898e-3       # Wien's displacement constant (m*K)

    # --- Given Measurements ---
    L_nm = 400.0  # Wavelength in nanometers
    T_K = 9000.0  # Temperature in Kelvin
    B_W = 1.2e15  # Spectral radiance in W/(m^2*sr*m)

    # Convert wavelength to meters for calculations
    L_m = L_nm * 1e-9

    # --- Analysis ---
    # For a given temperature, there is a maximum possible spectral radiance.
    # Let's calculate this peak radiance for T = 9000 K.
    
    # 1. Find the peak wavelength using Wien's displacement law: L_max = b / T
    L_max_m = b / T_K

    # 2. Calculate the peak spectral radiance B_peak at L_max_m and T_K
    exponent_peak = (h * c) / (L_max_m * k * T_K)
    B_peak = (2 * h * c**2) / (L_max_m**5 * (math.exp(exponent_peak) - 1))

    # 3. Compare the measured radiance B_W with the calculated peak radiance B_peak.
    # If B_W > B_peak, the measurements are physically inconsistent.
    if B_W > B_peak:
        # The measured radiance is higher than the maximum possible radiance for the given temperature.
        # This indicates a clear error. Given the astrophysical context for DB white dwarfs,
        # the temperature measurement is the most likely error.
        
        # We will now calculate the temperature T' that would produce the measured radiance B_W
        # at the given wavelength L_m by solving Planck's law for T.
        
        # Calculate the terms for the equation
        two_h_c_squared = 2 * h * c**2
        L_m_5 = L_m**5
        
        # Calculate the corrected temperature
        term_in_ln = two_h_c_squared / (B_W * L_m_5)
        T_corrected = (h * c) / (L_m * k * math.log(term_in_ln + 1))
        
        # As requested, output the equation with the numbers plugged in.
        print("The measured values are inconsistent. The measured radiance is physically impossible for the given temperature.")
        print("Calculating the corrected temperature T' assuming radiance and wavelength are correct.")
        print("Final Equation:")
        print(f"T' = ({h * c}) / ({L_m} * {k} * math.log(( {two_h_c_squared} ) / ( {B_W} * {L_m_5} ) + 1))")
        
        # Print the final answer in the required format Xv
        print(f"T{round(T_corrected)}")

    else:
        # If the measurement is physically possible, we check if it's reasonably close.
        exponent_b = (h * c) / (L_m * k * T_K)
        B_calc = (2 * h * c**2) / (L_m**5 * (math.exp(exponent_b) - 1))
        
        # Check if the relative difference is small (e.g., less than 10%)
        if abs(B_calc - B_W) / B_W < 0.1:
            print("0")
        else:
            # The values are inconsistent but not physically impossible.
            # Following the same logic, the temperature is the most likely error.
            term_in_ln = (2 * h * c**2) / (B_W * L_m**5)
            T_corrected = (h * c) / (L_m * k * math.log(term_in_ln + 1))
            print(f"T{round(T_corrected)}")

solve_star_measurement()