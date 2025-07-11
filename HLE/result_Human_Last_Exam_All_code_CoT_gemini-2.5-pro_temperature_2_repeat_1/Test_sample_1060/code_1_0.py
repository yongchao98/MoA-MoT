def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles available
    per engine revolution based on given parameters.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # MegaHertz

    # 1. Convert units to a common base (seconds)
    # ECU speed from MHz to Hertz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # Engine speed from RPM to RPS (revolutions per second)
    # This also represents the number of interrupts per second.
    revolutions_per_second = engine_rpm / 60

    # 2. Calculate the theoretical maximum number of cycles per revolution
    # This is the total cycles per second divided by revolutions per second.
    max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

    # 3. Print the final equation with the numbers used in the calculation
    print("The theoretical maximum number of cycles is calculated by dividing the ECU's cycles per second by the engine's revolutions per second.")
    print("Final Equation:")
    print(f"({int(ecu_speed_hz)} cycles/sec) / ({int(engine_rpm)} revs / 60 sec) = {int(max_cycles_per_revolution)} cycles")

# Run the calculation and print the result
calculate_ecu_cycles()