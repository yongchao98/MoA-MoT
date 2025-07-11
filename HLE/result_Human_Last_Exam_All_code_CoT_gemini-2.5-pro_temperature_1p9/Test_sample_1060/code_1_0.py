def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450   # MegaHertz

    # 1. Convert engine speed to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60

    # 2. Convert ECU speed to Hertz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # 3. Calculate the maximum number of cycles available per revolution.
    # This is the total cycles per second divided by the revolutions per second.
    max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

    # Print the explanation and the final equation with values
    print("To find the maximum cycles per revolution, we use the following equation:")
    print("Cycles per Revolution = (ECU Speed in Hz) / (Engine Speed in RPS)\n")
    print("Calculation:")
    print(f"{int(ecu_speed_hz):,} cycles/sec / {int(revolutions_per_second)} revs/sec = {int(max_cycles_per_revolution):,} cycles per revolution")


calculate_ecu_cycles()
<<<6000000>>>