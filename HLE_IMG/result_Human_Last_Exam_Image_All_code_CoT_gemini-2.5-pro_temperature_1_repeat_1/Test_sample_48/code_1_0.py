def estimate_martensite_start_temperature():
    """
    This script estimates the martensite start (Ms) temperature by analyzing
    the relationship between interpass temperature and residual stress in LTTE welds.
    """

    print("Step 1: Define the underlying principle of LTTE welds.")
    print("The volumetric expansion from the austenite-to-martensite phase transformation is used to generate compressive residual stress.")
    print("To be effective, this transformation must occur during the final cooling of the weld.")
    print("This requires the interpass temperature to be maintained ABOVE the martensite start (Ms) temperature.\n")

    # Data is qualitatively extracted from the 'Model' plots in Figure 1.
    # The color legend shows dark brown/black (-300 to -500 MPa) is high compression.
    # Yellow/light orange (0 to -200 MPa) is low/moderate compression.
    observations = {
        50: "Low compressive stress (light orange/yellow color in the weld center).",
        100: "High compressive stress (dark brown/black color, indicating values of -300 to -500 MPa).",
        150: "High compressive stress (dark brown/black color, similar to 100°C).",
        200: "Reduced compressive stress (stress is less compressive than at 100/150°C)."
    }

    print("Step 2: Analyze the observations from the stress maps.")
    for temp, desc in observations.items():
        print(f"- At an interpass temperature of {temp}°C: {desc}")

    print("\nStep 3: Determine the critical temperature range for Ms.")
    print("We are looking for the transition from ineffective to effective generation of compressive stress.")
    print("- The result at an interpass temperature of 50°C is suboptimal. This suggests that 50°C is likely below or too close to the Ms temperature.")
    print("- The result at an interpass temperature of 100°C shows highly effective generation of compressive stress. This suggests that 100°C is above the Ms temperature.")
    
    lower_bound = 50
    upper_bound = 100
    
    print(f"\nThe significant improvement in compressive stress occurs between {lower_bound}°C and {upper_bound}°C.")
    print(f"Therefore, the martensite start (Ms) temperature must lie within this range.\n")
    
    print(f"Conclusion: The estimated Ms temperature is between {lower_bound}°C and {upper_bound}°C.")
    print("This corresponds to answer choice E.")

# Run the analysis
estimate_martensite_start_temperature()