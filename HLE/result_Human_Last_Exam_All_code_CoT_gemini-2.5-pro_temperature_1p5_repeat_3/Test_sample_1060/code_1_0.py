def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450      # Megahertz

    # 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
    engine_rps = engine_rpm / 60

    # 2. Convert ECU speed from MHz to Hz (cycles per second)
    ecu_hz = ecu_mhz * 1_000_000

    # 3. Calculate the number of available cycles per revolution
    # This is the total cycles per second divided by the revolutions per second
    cycles_per_revolution = ecu_hz / engine_rps

    # Print the final equation showing all the numbers
    print(f"Engine speed: {engine_rpm} RPM")
    print(f"ECU processor speed: {ecu_mhz} MHz")
    print("\nCalculation:")
    print(f"First, we find the number of processor cycles per second: {ecu_mhz} MHz * 1,000,000 = {int(ecu_hz)} cycles/sec")
    print(f"Next, we find the engine revolutions per second: {engine_rpm} RPM / 60 = {int(engine_rps)} rev/sec")
    print("\nFinally, we divide cycles/sec by rev/sec to get cycles/rev:")
    print(f"{int(ecu_hz)} / {int(engine_rps)} = {int(cycles_per_revolution)}")
    print(f"\nThe theoretical maximum number of cycles is {int(cycles_per_revolution)}.")

calculate_ecu_cycles()
<<<6000000>>>