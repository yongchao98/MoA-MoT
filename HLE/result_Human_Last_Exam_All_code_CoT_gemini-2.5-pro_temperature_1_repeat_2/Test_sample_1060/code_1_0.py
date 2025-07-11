def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles available between interrupts.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_freq_mhz = 450   # ECU processor frequency in MegaHertz

    # Step 1: Calculate the time between interrupts
    # Convert RPM to Revolutions Per Second (RPS)
    rps = engine_rpm / 60.0
    # The time for one revolution is the time between interrupts
    time_per_interrupt_s = 1.0 / rps

    # Step 2: Calculate CPU cycles per second
    # Convert MHz to Hz (cycles per second)
    ecu_freq_hz = ecu_freq_mhz * 1_000_000

    # Step 3: Calculate the maximum number of cycles between interrupts
    max_cycles = ecu_freq_hz * time_per_interrupt_s

    # Output the explanation and the final calculation
    print("To find the maximum number of cycles, we first calculate the time available between interrupts.")
    print(f"Time per interrupt = 1 / ({engine_rpm} RPM / 60 s/min) = {time_per_interrupt_s:.6f} seconds")
    print("\nNext, we calculate how many cycles the ECU can perform per second.")
    print(f"ECU Speed = {ecu_freq_mhz} MHz * 1,000,000 = {int(ecu_freq_hz):,} cycles/second")
    print("\nFinally, we multiply the ECU speed by the time per interrupt to get the max cycles.")
    print("\n--- Final Equation ---")
    print(f"Max Cycles = {int(ecu_freq_hz):,} cycles/second * {time_per_interrupt_s:.6f} seconds")
    print(f"Max Cycles = {int(max_cycles):,}")

calculate_ecu_cycles()
<<<6000000>>>