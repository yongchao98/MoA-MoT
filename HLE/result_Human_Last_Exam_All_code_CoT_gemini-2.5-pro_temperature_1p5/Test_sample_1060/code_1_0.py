def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # Megahertz

    # 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60

    # 2. Convert ECU speed from MHz to cycles per second (Hz)
    cycles_per_second = ecu_speed_mhz * 1_000_000

    # 3. Calculate the maximum number of cycles per revolution
    # This is found by dividing the total cycles per second by the revolutions per second
    max_cycles_per_revolution = cycles_per_second / revolutions_per_second

    # Output the steps and the final calculation as requested
    print(f"Given an engine speed of {engine_rpm} RPM and an ECU speed of {ecu_speed_mhz} MHz.")
    print("-" * 30)
    print("Step 1: Calculate Revolutions Per Second (RPS)")
    print(f"{engine_rpm} RPM / 60 seconds/minute = {revolutions_per_second} RPS\n")

    print("Step 2: Calculate ECU Cycles Per Second (Hz)")
    print(f"{ecu_speed_mhz} MHz * 1,000,000 = {int(cycles_per_second)} Hz\n")

    print("Step 3: Calculate Max Cycles Per Revolution")
    print("This is the number of available CPU cycles between the once-per-revolution interrupts.")
    # The final equation with all numbers included
    print(f"Final Equation: {int(cycles_per_second)} cycles/sec / {int(revolutions_per_second)} revs/sec = {int(max_cycles_per_revolution)} cycles/rev\n")

    print(f"The theoretical maximum number of cycles you can run in the algorithm is: {int(max_cycles_per_revolution)}")

# Execute the function
calculate_ecu_cycles()
