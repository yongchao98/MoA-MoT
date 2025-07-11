def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles available per engine revolution.
    """
    # Step 1: Define the initial parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450      # Megahertz

    # Step 2: Convert units to be consistent (per second)
    # Convert ECU speed to Hertz (cycles per second)
    ecu_hz = ecu_mhz * 1_000_000

    # Convert engine speed to Revolutions Per Second (RPS)
    engine_rps = engine_rpm / 60

    # Step 3: Calculate the maximum cycles per revolution
    # The time window for calculation is one revolution.
    # Max cycles = Cycles per Second / Revolutions per Second = Cycles per Revolution
    max_cycles = ecu_hz / engine_rps

    # Print the explanation and the final equation
    print("Problem:")
    print(f"An engine is running at {engine_rpm} RPM with an ECU at {ecu_mhz} MHz.")
    print("An interrupt is triggered every revolution. Find the max CPU cycles per interrupt.\n")

    print("Step 1: Convert Engine RPM to Revolutions Per Second (RPS)")
    print(f"{engine_rpm} RPM / 60 seconds/minute = {int(engine_rps)} RPS\n")

    print("Step 2: Convert ECU MHz to Hertz (cycles per second)")
    print(f"{ecu_mhz} MHz * 1,000,000 = {int(ecu_hz)} Hz\n")

    print("Step 3: Calculate the maximum cycles available per revolution")
    print("This is found by dividing the ECU's cycles per second by the engine's revolutions per second.\n")

    print("Final Equation:")
    print("ECU Speed (Hz) / Engine Speed (RPS) = Max Cycles per Revolution")
    print(f"{int(ecu_hz)} / {int(engine_rps)} = {int(max_cycles)}")


calculate_ecu_cycles()