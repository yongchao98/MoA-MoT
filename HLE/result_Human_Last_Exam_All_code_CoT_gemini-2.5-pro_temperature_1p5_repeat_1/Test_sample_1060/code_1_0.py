def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles for an ECU algorithm
    based on engine speed and processor frequency.
    """
    # Given parameters
    engine_rpm = 4500  # Engine speed in Revolutions Per Minute
    ecu_mhz = 450      # ECU processor speed in MegaHertz

    # --- Step 1 & 2: Calculate time per revolution ---
    # Convert RPM to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60.0
    # The time for one revolution is the inverse of RPS
    seconds_per_revolution = 1.0 / revolutions_per_second

    # --- Step 3: Convert ECU speed to Hz ---
    # 1 MHz = 1,000,000 Hz (cycles per second)
    ecu_hz = ecu_mhz * 1_000_000

    # --- Step 4: Calculate maximum cycles per revolution ---
    # Cycles per revolution = (Cycles / Second) * (Seconds / Revolution)
    max_cycles = ecu_hz * seconds_per_revolution

    # --- Output the results ---
    print("This script calculates the maximum CPU cycles available per engine revolution.")
    print(f"Engine Speed: {engine_rpm} RPM")
    print(f"ECU Processor Speed: {ecu_mhz} MHz\n")

    # The user requested to output each number in the final equation.
    # The equation is: (ECU Speed in Hz) / (Revolutions per Second)
    # Revolutions per Second = engine_rpm / 60
    print("The final calculation is:")
    print(f"({int(ecu_hz):,}) cycles/sec / ({int(engine_rpm)} revs / 60 sec) = {int(max_cycles):,} cycles per revolution")

if __name__ == "__main__":
    calculate_ecu_cycles()
