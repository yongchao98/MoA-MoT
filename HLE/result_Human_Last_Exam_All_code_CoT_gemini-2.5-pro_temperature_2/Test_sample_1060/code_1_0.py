def calculate_cycles_per_revolution():
    """
    Calculates the theoretical maximum ECU cycles per engine revolution.
    """
    # --- Given Parameters ---
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450    # Processor speed in MegaHertz

    # --- Step 1: Convert Engine Speed to Revolutions Per Second (RPS) ---
    seconds_per_minute = 60
    engine_rps = engine_rpm / seconds_per_minute

    # --- Step 2: Convert ECU Speed to Cycles Per Second (Hz) ---
    mhz_to_hz_multiplier = 1_000_000
    ecu_speed_hz = ecu_speed_mhz * mhz_to_hz_multiplier

    # --- Step 3: Calculate Maximum Cycles per Revolution ---
    # This is the total number of cycles per second divided by the revolutions per second.
    max_cycles = ecu_speed_hz / engine_rps

    # --- Output the final equation with numbers ---
    print(f"To find the cycles per revolution, we divide the processor's cycles per second by the engine's revolutions per second.")
    print(f"Processor Speed in Hz: {ecu_speed_mhz} MHz * {mhz_to_hz_multiplier} = {int(ecu_speed_hz)} cycles/sec")
    print(f"Engine Speed in RPS: {engine_rpm} RPM / {seconds_per_minute} sec/min = {int(engine_rps)} revs/sec")
    print(f"\nFinal Equation:")
    print(f"{int(ecu_speed_hz)} / {int(engine_rps)} = {int(max_cycles)} cycles")

if __name__ == "__main__":
    calculate_cycles_per_revolution()