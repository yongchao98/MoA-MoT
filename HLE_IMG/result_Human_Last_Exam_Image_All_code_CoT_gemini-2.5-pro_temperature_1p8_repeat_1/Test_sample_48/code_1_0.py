def estimate_martensite_start_temperature():
    """
    This function analyzes the relationship between interpass temperature and residual
    stress in LTTE welds to estimate the Martensite Start (Ms) temperature.
    """
    
    print("Step 1: Understand the principle of LTTE welds.")
    print("Low Transformation Temperature Expansion (LTTE) welds create compressive residual stress.")
    print("This happens because of a volume expansion when the material transforms to martensite at a low temperature (Ms).\n")
    
    print("Step 2: Analyze the effect of Interpass Temperature (T_interpass).")
    print(" - If T_interpass > Ms: The transformation occurs only after all welding is done. This creates high compressive stress.")
    print(" - If T_interpass < Ms: The transformation occurs between weld passes. The heat from subsequent passes removes the beneficial compressive stress.\n")

    print("Step 3: Observe the stress maps from the provided image.")
    print(" - For T_interpass = 200 C, 150 C, and 100 C:")
    print("   The weld shows high compressive stresses (dark-colored regions, -300 to -500 MPa).")
    print("   This implies that the material was kept above its Ms temperature between passes.")
    print("   Therefore, the Martensite start temperature (Ms) must be less than 100 C.\n")

    print(" - For T_interpass = 50 C:")
    print("   The high compressive stresses have disappeared. The stress is near zero (white/light-yellow region).")
    print("   This implies the beneficial effect was lost because the material cooled below Ms between passes.")
    print("   Therefore, the Martensite start temperature (Ms) must be greater than 50 C.\n")

    print("Step 4: Conclude the estimated range for Ms.")
    low_bound = 50
    high_bound = 100
    print(f"Combining the observations, the Ms temperature must be greater than {low_bound} C and less than {high_bound} C.")
    print(f"Final Estimated Range for Ms: {low_bound}°C - {high_bound}°C")

if __name__ == "__main__":
    estimate_martensite_start_temperature()