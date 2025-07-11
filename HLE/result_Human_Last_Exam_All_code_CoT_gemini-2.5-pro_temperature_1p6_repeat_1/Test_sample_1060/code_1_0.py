def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450      # Megahertz

    # Step 1: Convert RPM to Revolutions Per Second (RPS)
    seconds_in_minute = 60
    engine_rps = engine_rpm / seconds_in_minute

    # Step 2: Convert ECU speed from MHz to Hz (cycles per second)
    mhz_to_hz_conversion = 1_000_000
    ecu_hz = ecu_mhz * mhz_to_hz_conversion

    # Step 3: Calculate the maximum cycles per revolution (interrupt)
    # This is done by dividing the total cycles per second by the revolutions per second.
    cycles_per_interrupt = ecu_hz / engine_rps

    # Output the explanation and the final equation
    print(f"The engine runs at {engine_rpm} RPM and the ECU runs at {ecu_mhz} MHz.")
    print("An interrupt occurs once per revolution.")
    print("To find the cycles available per interrupt, we divide the ECU's cycles per second by the engine's revolutions per second.\n")

    print("Final Calculation:")
    # Print the full equation with each number clearly formatted
    print(f"({ecu_mhz} MHz * {mhz_to_hz_conversion:,}) cycles/sec / ({engine_rpm} RPM / {seconds_in_minute} sec/min) revs/sec")
    print(f"= {int(ecu_hz):,} cycles/sec / {int(engine_rps)} revs/sec")
    print(f"= {int(cycles_per_interrupt):,} cycles per interrupt")

calculate_ecu_cycles()