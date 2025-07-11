import math

def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles for an ECU algorithm.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_mhz = 450      # Megahertz

    # Step 1: Convert engine speed to Revolutions Per Second (RPS)
    engine_rps = engine_rpm / 60

    # Step 2: Convert ECU speed to cycles per second (Hz)
    ecu_hz = ecu_mhz * 1_000_000

    # Step 3: Calculate the maximum cycles per revolution
    # This is the theoretical budget available between interrupts
    max_cycles = ecu_hz / engine_rps

    # Ensure the result is a clean integer
    final_cycles_int = int(max_cycles)
    
    # Print the full equation with each number, as requested
    print("The final equation is:")
    # Format numbers with commas for readability
    print(f"{int(ecu_hz):,} cycles/sec / {int(engine_rps)} revs/sec = {final_cycles_int:,} cycles/rev")
    
    # Output the final answer in the specified format
    print(f"\n<<<{final_cycles_int}>>>")

# Execute the function
calculate_ecu_cycles()