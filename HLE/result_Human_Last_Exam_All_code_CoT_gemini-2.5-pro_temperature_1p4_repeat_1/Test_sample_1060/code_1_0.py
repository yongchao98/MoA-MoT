def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles available
    between engine revolutions for an ECU.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # Megahertz

    # 1. Convert engine RPM to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60.0

    # 2. Convert ECU speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # 3. Calculate the maximum number of cycles per revolution (i.e., per interrupt)
    # This is done by dividing the total cycles per second by the revolutions per second.
    max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

    # 4. Print the final equation and result
    print(f"The theoretical maximum number of cycles is calculated as follows:")
    print(f"ECU Speed (Hz) / Revolutions per Second (RPS)")
    print(f"Which is: {int(ecu_speed_hz):,} / ({int(engine_rpm)} / 60)")
    print(f"Result: {int(max_cycles_per_revolution):,} cycles")

calculate_ecu_cycles()
<<<6000000>>>