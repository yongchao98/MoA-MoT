def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles for an ECU algorithm
    based on engine RPM and processor speed.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450       # ECU speed in MegaHertz

    # 1. Convert engine RPM to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60

    # 2. Convert ECU speed from MHz to Hz (cycles per second)
    cycles_per_second = ecu_mhz * 1_000_000

    # 3. Calculate the maximum number of cycles available per revolution.
    # This is the total cycles per second divided by the revolutions per second.
    max_cycles_per_revolution = cycles_per_second / revolutions_per_second

    # Print the final equation and the result
    # The equation is: (CPU Cycles per Second) / (Revolutions per Second)
    print(f"The equation to find the maximum cycles per revolution is:")
    print(f"({int(cycles_per_second)}) / ({int(engine_rpm)} / 60) = {int(max_cycles_per_revolution)}")
    print(f"\nTherefore, the theoretical maximum number of cycles you can run in the algorithm is {int(max_cycles_per_revolution)}.")

calculate_ecu_cycles()