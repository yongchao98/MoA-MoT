def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles for a spark plug
    algorithm given engine RPM and ECU speed.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450  # MegaHertz

    # Convert units
    # Processor speed in cycles per second (Hz)
    ecu_hz = ecu_mhz * 1000000
    # Engine speed in revolutions per second (RPS)
    engine_rps = engine_rpm / 60

    # Calculate the maximum cycles available per revolution
    # This is found by dividing the cycles per second by revolutions per second
    if engine_rps > 0:
        max_cycles = ecu_hz / engine_rps
    else:
        max_cycles = 0

    # Output the full equation and the result
    print(f"Given an engine speed of {engine_rpm} RPM and an ECU speed of {ecu_mhz} MHz:")
    print("The maximum number of cycles per revolution is calculated as:")
    print(f"({ecu_hz:,} cycles/sec) / ({engine_rpm} rev/min / 60 sec/min) = {int(max_cycles):,} cycles")

calculate_ecu_cycles()

# The final numerical answer
print("\n<<<6000000>>>")