def calculate_ecu_cycles():
    """
    Calculates the maximum number of ECU cycles available per engine revolution.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450      # ECU processor speed in MegaHertz

    # Step 1: Calculate time per revolution in seconds
    # Time per Revolution = 60 seconds / RPM
    seconds_per_revolution = 60 / engine_rpm

    # Step 2: Calculate ECU cycles per second (Hertz)
    # 1 MHz = 1,000,000 Hz
    ecu_hz = ecu_mhz * 1_000_000

    # Step 3: Calculate the maximum cycles available per revolution
    # This is the number of cycles per second multiplied by the seconds per revolution
    max_cycles = ecu_hz * seconds_per_revolution

    # Output the explanation and the final equation
    print("This script calculates the maximum number of ECU cycles per engine revolution based on the given data.\n")
    print(f"1. Engine speed is {engine_rpm} RPM.")
    print(f"2. ECU processor runs at {ecu_mhz} MHz ({int(ecu_hz)} cycles per second).")
    print(f"3. An interrupt occurs every revolution, giving the algorithm a time window of {seconds_per_revolution:.4f} seconds.\n")
    
    print("The final calculation is:")
    print(f"({ecu_mhz} * 1000000 cycles/sec) * (60 sec / {engine_rpm} rev) = {int(max_cycles)} cycles/rev")


calculate_ecu_cycles()
<<<6000000>>>