def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of ECU cycles per engine revolution.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_frequency_mhz = 450  # Megahertz

    # Step 1: Convert RPM to Revolutions Per Second (RPS)
    revolutions_per_second = engine_rpm / 60.0

    # Step 2: Convert ECU frequency from MHz to Hz (cycles per second)
    ecu_frequency_hz = ecu_frequency_mhz * 1_000_000

    # Step 3: Calculate the maximum number of cycles per revolution (per interrupt)
    # This is done by dividing the total cycles per second by the revolutions per second.
    max_cycles_per_revolution = ecu_frequency_hz / revolutions_per_second

    # Output the result, showing each number in the equation
    print("This script calculates the maximum ECU cycles available per engine revolution.")
    print("-" * 60)
    print(f"Engine Speed: {engine_rpm} RPM")
    print(f"ECU Processor Speed: {ecu_frequency_mhz} MHz\n")

    print("Calculation:")
    print(f"({ecu_frequency_hz:,} cycles/sec) / ({engine_rpm} rev/min / 60 sec/min) = {int(max_cycles_per_revolution):,} cycles/rev")

    # Print a more direct version of the final equation
    print(f"\nFinal Equation using values in Hz and RPS:")
    print(f"{int(ecu_frequency_hz):,} / {int(revolutions_per_second):,} = {int(max_cycles_per_revolution):,}")
    print("-" * 60)
    print(f"The theoretical maximum number of cycles per revolution is {int(max_cycles_per_revolution):,}.")

calculate_ecu_cycles()

# The final numerical answer
print(f"\n<<<{int(6000000)}>>>")