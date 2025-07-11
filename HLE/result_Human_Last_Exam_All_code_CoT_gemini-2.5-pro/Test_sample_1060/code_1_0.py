def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450       # MegaHertz

    # 1. Convert engine RPM to Revolutions Per Second (RPS)
    seconds_in_minute = 60
    engine_rps = engine_rpm / seconds_in_minute

    # 2. Convert ECU MHz to Hz (cycles per second)
    hz_in_mhz = 1_000_000
    ecu_hz = ecu_mhz * hz_in_mhz

    # 3. Calculate the number of ECU cycles per revolution
    # This is the theoretical maximum number of cycles available for the algorithm
    max_cycles_per_revolution = ecu_hz / engine_rps

    # Output the final equation with all the numbers
    print(f"Engine speed: {engine_rpm} RPM")
    print(f"ECU processor speed: {ecu_mhz} MHz")
    print("\nTo find the maximum cycles per revolution, we use the following calculation:")
    print(f"({ecu_mhz} * {hz_in_mhz}) cycles/sec / ({engine_rpm} rev/min / {seconds_in_minute} sec/min) = {int(max_cycles_per_revolution)} cycles/rev")

calculate_ecu_cycles()
<<<6000000>>>