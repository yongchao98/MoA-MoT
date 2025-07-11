def calculate_cycles_per_revolution():
    """
    Calculates the theoretical maximum CPU cycles available per engine revolution.
    """
    # Step 1: Establish given values
    engine_rpm = 4500
    ecu_speed_mhz = 450
    
    # Step 2: Standardize units for calculation
    # Convert ECU speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000
    # Define seconds in a minute for RPM to RPS conversion
    seconds_per_minute = 60

    # Step 3: Calculate the cycles per revolution
    # Engine speed in Revolutions Per Second (RPS)
    engine_rps = engine_rpm / seconds_per_minute
    # Total cycles per revolution = (total cycles per second) / (revolutions per second)
    cycles_per_revolution = ecu_speed_hz / engine_rps

    # Step 4: Display the result, including the full equation as requested
    print(f"The theoretical maximum number of cycles available to the ECU algorithm per revolution is calculated as follows:")
    print(f"({int(ecu_speed_hz)} cycles/sec) / ({engine_rpm} rev/min / {seconds_per_minute} sec/min) = {int(cycles_per_revolution)} cycles")

calculate_cycles_per_revolution()

<<<6000000>>>