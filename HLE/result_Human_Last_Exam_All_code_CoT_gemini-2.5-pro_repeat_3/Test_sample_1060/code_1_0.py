def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles available per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450 # MegaHertz

    # Conversion factors
    seconds_per_minute = 60
    hertz_per_megahertz = 1_000_000

    # Step 1: Convert engine speed to Revolutions Per Second (RPS)
    engine_rps = engine_rpm / seconds_per_minute

    # Step 2: Convert ECU speed to Hertz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * hertz_per_megahertz

    # Step 3: Calculate the maximum available cycles per revolution
    # (cycles/second) / (revolutions/second) = cycles/revolution
    max_cycles = ecu_speed_hz / engine_rps

    print(f"An engine running at {engine_rpm} RPM completes {engine_rps} revolutions per second.")
    print(f"An ECU running at {ecu_speed_mhz} MHz can perform {int(ecu_speed_hz):,} cycles per second.")
    print("\nTo find the maximum cycles per revolution, we divide the ECU's cycles per second by the engine's revolutions per second:")
    
    # Print the final equation with all the numbers
    print(f"\n{int(ecu_speed_hz):,} / {int(engine_rps)} = {int(max_cycles):,}")
    
    # The result is the theoretical maximum number of cycles for the algorithm
    print(f"\nTherefore, the theoretical maximum number of cycles you can run in the algorithm is {int(max_cycles):,}.")

if __name__ == "__main__":
    calculate_ecu_cycles()