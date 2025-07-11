def estimate_martensite_start_temperature():
    """
    Analyzes the provided stress maps to estimate the martensite start temperature (Ms).
    """

    print("Step 1: Understand the principle of LTTE fillers.")
    print("Low Transformation Temperature Expansion (LTTE) fillers create compressive residual stress.")
    print("This happens because of a volume expansion during the austenite-to-martensite phase transformation, which counteracts thermal shrinkage during cooling.")
    print("This beneficial effect only occurs if the weld material cools below its Martensite start (Ms) temperature.\n")

    print("Step 2: Analyze the residual stress maps at different interpass temperatures.")
    print("The experiment uses interpass temperature to control the minimum temperature the weld reaches between passes.")

    temp_100 = 100
    print(f" - At an interpass temperature of {temp_100}°C (and 50°C), the maps show high compressive stress in the weld.")
    print("   This indicates that the weld cooled below the Ms temperature, allowing the transformation to occur and generate compression.\n")

    temp_150 = 150
    print(f" - At an interpass temperature of {temp_150}°C, the compressive stress is significantly reduced.")
    print("   This is the key observation. It means that keeping the weld from cooling below 150°C prevents the transformation.")
    print(f"   Therefore, the Ms temperature must be below {temp_150}°C.\n")
    
    temp_200 = 200
    print(f" - At an interpass temperature of {temp_200}°C, the weld shows tensile stress, confirming the Ms temperature is well below {temp_200}°C.\n")

    print("Step 3: Conclude the estimated range for Ms.")
    print(f"The dramatic change in residual stress behavior occurs between an interpass temperature of {temp_100}°C and {temp_150}°C.")
    print(f"This strongly suggests that the Martensite start (Ms) temperature of the filler material is in the range of {temp_100}°C to {temp_150}°C.")

estimate_martensite_start_temperature()