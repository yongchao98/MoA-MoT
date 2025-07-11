def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450 # MegaHertz

    # 1. Convert ECU speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # 2. Convert engine RPM to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60

    # 3. Calculate the number of available cycles per revolution
    cycles_per_revolution = ecu_speed_hz / revolutions_per_second

    # Print the equation and the final result
    print("ECU Speed (Hz) / Engine Speed (RPS) = Max Cycles per Revolution")
    print(f"({int(ecu_speed_hz):,d} cycles/sec) / ({int(engine_rpm):,d} RPM / 60 sec/min) = {int(cycles_per_revolution):,d} cycles")

calculate_ecu_cycles()