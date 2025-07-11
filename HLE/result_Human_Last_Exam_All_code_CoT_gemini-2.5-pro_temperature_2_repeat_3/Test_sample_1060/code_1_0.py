def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450       # MegaHertz

    # 1. Convert engine RPM to Revolutions Per Second (RPS)
    rpm_to_rps_divisor = 60
    engine_rps = engine_rpm / rpm_to_rps_divisor

    # 2. Convert ECU MHz to Hz (cycles per second)
    mhz_to_hz_multiplier = 1_000_000
    ecu_hz = ecu_mhz * mhz_to_hz_multiplier

    # 3. Calculate the maximum number of cycles per revolution
    # This is found by dividing the processor's cycles per second by the engine's revolutions per second.
    max_cycles_per_revolution = ecu_hz / engine_rps

    # Print the explanation and the final equation with the numbers
    print(f"The engine runs at {engine_rpm} RPM, which is {engine_rps} revolutions per second.")
    print(f"The ECU runs at {ecu_mhz} MHz, which is {ecu_hz:,} cycles per second.")
    print("The theoretical maximum number of cycles per revolution is calculated as:")
    print(f"({ecu_hz:,} cycles/sec) / ({engine_rpm} rev/min / {rpm_to_rps_divisor} sec/min) = {int(max_cycles_per_revolution):,} cycles")

calculate_ecu_cycles()
<<<6000000>>>