def calculate_ecu_cycles():
    """
    Calculates the maximum number of CPU cycles available for an algorithm
    between engine revolutions.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # Megahertz

    # Step 1: Convert ECU speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # Step 2: Convert engine speed from RPM to RPS (revolutions per second)
    # The number of interrupts per second is equal to the RPS
    interrupts_per_second = engine_rpm / 60

    # Step 3: Calculate the maximum cycles available between interrupts (per revolution)
    # This is the processor speed (cycles/sec) divided by the number of interrupts (interrupts/sec)
    if interrupts_per_second == 0:
        max_cycles = float('inf')
    else:
        max_cycles = ecu_speed_hz / interrupts_per_second

    # --- Output the result ---
    # We want to display the calculation clearly.
    # The equation is: (Processor Speed in Hz) / (Engine Speed in RPS)
    # Engine Speed in RPS is calculated as (Engine RPM / 60)
    print("ECU Processor Speed: {:,} MHz = {:,} cycles/second".format(ecu_speed_mhz, ecu_speed_hz))
    print("Engine Speed: {:,} RPM".format(engine_rpm))
    print("\nAn interrupt fires every revolution.")
    print("To find the maximum cycles per revolution, we use the formula:")
    print("Max Cycles = (Processor Speed in Hz) / (Revolutions per Second)\n")

    # Print the final equation with all numbers
    print("Calculation:")
    print("({:,} cycles/sec) / (({:,} rev/min) / 60 sec/min) = {:,} cycles/revolution".format(
        ecu_speed_hz,
        engine_rpm,
        int(max_cycles)
    ))

calculate_ecu_cycles()
<<<6000000>>>