def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles for an ECU
    between interrupts that occur once per engine revolution.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    processor_speed_mhz = 450  # Megahertz

    # --- Step 1: Standardize Units ---
    # Convert processor speed from MHz to Hz (cycles per second)
    processor_speed_hz = processor_speed_mhz * 1_000_000

    # Convert engine speed from RPM to RPS (revolutions per second)
    engine_rps = engine_rpm / 60

    # --- Step 2: Calculate Cycles Per Revolution ---
    # Divide total cycles per second by revolutions per second
    max_cycles = processor_speed_hz / engine_rps

    # --- Step 3: Output the results ---
    print("This calculation determines the number of processor cycles available for every single revolution of the engine.")
    print("\nFirst, we convert the given values to a common time unit of seconds:")
    print(f"Processor Speed: {processor_speed_mhz} MHz = {int(processor_speed_hz):,} cycles/second")
    print(f"Engine Speed: {engine_rpm} RPM = {engine_rps} revolutions/second")
    
    print("\nThe final equation is the processor speed (cycles/sec) divided by the engine speed (revs/sec):")
    print(f"({processor_speed_mhz} * 1,000,000) / ({engine_rpm} / 60)")
    print(f"= {int(processor_speed_hz):,} / {int(engine_rps)}")
    
    print(f"\nResult:")
    print(f"The theoretical maximum number of cycles available is {int(max_cycles):,}")

# Run the calculation
calculate_ecu_cycles()
