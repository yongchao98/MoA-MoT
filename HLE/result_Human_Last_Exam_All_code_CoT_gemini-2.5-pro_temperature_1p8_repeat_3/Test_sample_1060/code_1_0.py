def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles for an ECU
    algorithm between interrupts tied to engine revolutions.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # MegaHertz

    # 1. Convert ECU speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # 2. Calculate the time for one revolution in seconds
    # time = 60 seconds / revolutions_per_minute
    time_per_revolution_s = 60 / engine_rpm

    # 3. Calculate the total cycles available in that time window
    # cycles = cycles_per_second * seconds
    max_cycles = ecu_speed_hz * time_per_revolution_s

    print(f"Engine speed: {engine_rpm} RPM")
    print(f"ECU processor speed: {ecu_speed_mhz} MHz ({int(ecu_speed_hz):,} cycles/sec)")
    print("-" * 40)
    print(f"Time per revolution: 60 seconds / {engine_rpm} RPM = {time_per_revolution_s:.4f} seconds")
    print(f"Calculation for maximum cycles per revolution:")
    print(f"Final Equation: {int(ecu_speed_hz)} cycles/sec * {time_per_revolution_s:.4f} sec = {int(max_cycles):,} cycles")
    print("-" * 40)
    print(f"The theoretical maximum number of cycles is {int(max_cycles):,}.")

calculate_ecu_cycles()
<<<6000000>>>