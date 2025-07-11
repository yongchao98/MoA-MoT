import math

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's magnitude to a given precision.
    """
    # Explain the plan and assumptions
    print("This script calculates the required exposure time for a given photometric precision.")
    print("The calculation proceeds in the following steps:")
    print("1. Determine the required Signal-to-Noise Ratio (SNR) from the magnitude uncertainty.")
    print("2. Calculate the total number of photons needed to achieve this SNR.")
    print("3. Find the photon flux from the target star based on its magnitude.")
    print("4. Compute the photon collection rate for the given telescope.")
    print("5. Calculate the final exposure time by dividing the total photons by the collection rate.\n")

    # Define the given parameters and constants
    m_star = 20.0
    delta_m = 0.01
    D = 1.0  # Telescope diameter in meters

    # Define physical and astronomical constants/assumptions
    # The B-band zero-point flux is based on a spectral flux of 1500 photons/s/cm²/Å
    # and a filter bandwidth of 940 Å (94 nm). This is a standard value for a star
    # with B-band magnitude m_B = 0, observed from above the atmosphere.
    phi_0 = 1500 * 940 * 1e4  # photons/s/m^2
    
    # For this idealized problem, we assume 100% total system efficiency.
    efficiency = 1.0

    print("--- Input Parameters and Assumptions ---")
    print(f"Star B-band Magnitude (M_B): {m_star}")
    print(f"Desired Magnitude Accuracy (delta_m): +/- {delta_m}")
    print(f"Telescope Diameter (D): {D} m")
    print(f"Assumed B-band Zero-Point (m_B=0) Flux: {phi_0:.2e} photons/s/m^2")
    print(f"Assumed Overall Efficiency: {efficiency*100}%\n")

    # Step 1: Calculate the required Signal-to-Noise Ratio (SNR)
    # The relation is delta_m ≈ (2.5 / ln(10)) / SNR
    snr = (2.5 / math.log(10)) / delta_m
    print("--- Calculation Steps ---")
    print(f"1. Required Signal-to-Noise Ratio (SNR) = (2.5 / ln(10)) / {delta_m} = {snr:.2f}")

    # Step 2: Calculate the total number of photons (N) needed
    # For photon shot noise, SNR = sqrt(N)
    N = snr**2
    print(f"2. Total Photons Needed (N) = SNR^2 = {snr:.2f}^2 = {N:.2f}")

    # Step 3: Calculate the photon flux from the star
    # Flux = Phi_0 * 10^(-0.4 * m)
    phi_star = phi_0 * (10**(-0.4 * m_star))
    print(f"3. Photon Flux from Star = {phi_0:.2e} * 10^(-0.4 * {m_star:.1f}) = {phi_star:.2f} photons/s/m^2")

    # Step 4: Calculate the photon collection rate (R)
    # Rate = Flux * Telescope Area * Efficiency
    area = math.pi * (D / 2)**2
    rate = phi_star * area * efficiency
    print(f"4. Telescope Collection Rate (R) = Flux * Area = {phi_star:.2f} * {area:.4f} = {rate:.2f} photons/s")

    # Step 5: Calculate the final exposure time (t)
    time_s = N / rate
    print("\n--- Final Exposure Time Calculation ---")
    print(f"Exposure Time = Total Photons / Collection Rate")
    # Output the final equation with numbers
    print(f"              = {N:.0f} photons / {rate:.2f} photons/s")
    print(f"              = {time_s:.2f} seconds")

    # Round the final answer to the nearest integer
    time_rounded = round(time_s)
    print(f"\nThe required exposure time, rounded to the nearest integer, is {time_rounded} seconds.")

    return time_rounded

# Run the calculation and print the final answer in the required format
final_answer = calculate_exposure_time()
print(f"\n<<<{final_answer}>>>")
