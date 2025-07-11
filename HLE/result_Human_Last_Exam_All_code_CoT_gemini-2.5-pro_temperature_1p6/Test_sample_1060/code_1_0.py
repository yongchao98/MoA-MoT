def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles available
    for an algorithm between interrupts.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450   # ECU processor speed in MegaHertz

    # 1. Convert engine RPM to Revolutions Per Second (RPS)
    # This is the interrupt frequency since there's one interrupt per revolution.
    revolutions_per_second = engine_rpm / 60

    # 2. Convert ECU speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # 3. Calculate the maximum number of cycles available per revolution (between interrupts).
    # This is (cycles per second) / (revolutions per second).
    max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

    # Print the final equation with the numbers and the result
    print("The final equation to solve for the maximum available CPU cycles per revolution is:")
    # Using int() for cleaner display of whole numbers in the equation
    print(f"{int(ecu_speed_hz)} / ({int(engine_rpm)} / 60) = {int(max_cycles_per_revolution)}")

calculate_ecu_cycles()
<<<6000000>>>