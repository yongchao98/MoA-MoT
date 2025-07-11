def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles per engine revolution.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # Megahertz

    # --- Step 1: Convert Engine RPM to Revolutions Per Second (RPS) ---
    # This tells us how many interrupts occur per second.
    revolutions_per_second = engine_rpm / 60

    # --- Step 2: Convert ECU Speed from MHz to Hz (cycles per second) ---
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # --- Step 3: Calculate the maximum cycles per interrupt (revolution) ---
    # This is the total cycles per second divided by the interrupts per second.
    max_cycles_per_interrupt = ecu_speed_hz / revolutions_per_second

    # --- Output the results ---
    print(f"Engine Speed: {engine_rpm} RPM")
    print(f"ECU Processor Speed: {ecu_speed_mhz} MHz")
    print("-" * 30)
    print("Calculation:")
    print(f"({ecu_speed_mhz} * 1000000) cycles/sec / ({engine_rpm} / 60) revolutions/sec")
    print("=")
    print(f"{int(ecu_speed_hz)} cycles/sec / {int(revolutions_per_second)} revolutions/sec")
    print("=")
    print(f"{int(max_cycles_per_interrupt)} cycles/revolution")
    print("-" * 30)
    print("The theoretical maximum number of cycles available between each interrupt is:")
    print(int(max_cycles_per_interrupt))

if __name__ == "__main__":
    calculate_ecu_cycles()