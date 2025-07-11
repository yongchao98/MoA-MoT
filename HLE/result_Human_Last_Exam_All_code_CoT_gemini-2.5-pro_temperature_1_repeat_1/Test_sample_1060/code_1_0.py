def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # revolutions per minute
    ecu_speed_mhz = 450  # megahertz

    # --- Step 1: Convert Engine RPM to Revolutions per Second (RPS) ---
    seconds_per_minute = 60
    revolutions_per_second = engine_rpm / seconds_per_minute

    # --- Step 2: Convert ECU Speed to Hertz (Hz) ---
    mhz_to_hz_multiplier = 1_000_000
    ecu_speed_hz = ecu_speed_mhz * mhz_to_hz_multiplier

    # --- Step 3: Calculate Cycles per Revolution ---
    # This is the number of cycles the CPU can execute in the time it takes for one revolution.
    max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

    # --- Output the final equation and result ---
    print(
        "The theoretical maximum number of cycles is calculated by dividing the ECU speed in Hz by the engine speed in revolutions per second."
    )
    print(
        f"Equation: ({int(ecu_speed_hz)} cycles/sec) / ({int(engine_rpm)} rev/min / {seconds_per_minute} sec/min)"
    )
    print(f"Result: {int(max_cycles_per_revolution):,} cycles per revolution")

calculate_ecu_cycles()
<<<6000000>>>