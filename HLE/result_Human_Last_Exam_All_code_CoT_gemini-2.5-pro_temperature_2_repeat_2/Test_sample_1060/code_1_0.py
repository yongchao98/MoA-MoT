def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    rpm = 4500  # Engine speed in Revolutions Per Minute
    ecu_speed_mhz = 450  # ECU processor speed in MegaHertz

    # 1. Convert processor speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # 2. Convert engine speed from RPM to RPS (revolutions per second)
    revolutions_per_second = rpm / 60

    # 3. Calculate the maximum number of cycles per revolution
    # This is done by dividing the total cycles per second by the revolutions per second.
    max_cycles = ecu_speed_hz / revolutions_per_second
    
    # Ensure the final result is a whole number
    max_cycles = int(max_cycles)

    print("To find the maximum cycles per revolution, we divide the ECU's cycles per second by the engine's revolutions per second.")
    print("\nFinal Equation:")
    print(f"({ecu_speed_mhz} * 1000000) cycles/sec / ({rpm} / 60) revs/sec = {max_cycles} cycles/rev")

calculate_ecu_cycles()