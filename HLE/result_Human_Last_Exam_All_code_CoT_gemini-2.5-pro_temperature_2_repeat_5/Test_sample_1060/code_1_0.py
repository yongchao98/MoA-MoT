def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles an ECU can execute
    for its spark plug algorithm between interrupts.
    """
    # Given parameters
    engine_rpm = 4500  # Revolutions Per Minute
    cpu_speed_mhz = 450   # Processor speed in MegaHertz

    # Step 1: Convert engine speed from RPM to Revolutions Per Second (RPS)
    # There are 60 seconds in a minute.
    # RPS = RPM / 60
    revolutions_per_second = engine_rpm / 60

    # Step 2: Convert CPU speed from MHz to Hz (cycles per second)
    # 1 MHz = 1,000,000 Hz
    # Cycles per Second = MHz * 1,000,000
    cpu_cycles_per_second = cpu_speed_mhz * 1_000_000

    # Step 3: Calculate the maximum cycles available per revolution.
    # An interrupt occurs once per revolution, giving the ECU a time window
    # to perform its calculations.
    # Cycles per Revolution = (Cycles per Second) / (Revolutions per Second)
    max_cycles_per_revolution = cpu_cycles_per_second / revolutions_per_second

    # Output the explanation and the final equation with the numbers
    print("This script calculates the theoretical maximum CPU cycles available for an algorithm between engine revolutions.")
    print(f"Given an engine speed of {engine_rpm} RPM and an ECU speed of {cpu_speed_mhz} MHz.\n")
    print("The final calculation is the number of CPU cycles per second divided by the number of revolutions per second.")
    
    # Print the equation with the actual numbers used
    print("\nFinal Equation:")
    print(f"{int(cpu_cycles_per_second)} / {int(revolutions_per_second)} = {int(max_cycles_per_revolution)}")
    
    print(f"\nThis means the theoretical maximum number of cycles the ECU can run for the algorithm per revolution is {int(max_cycles_per_revolution):,}.")

if __name__ == "__main__":
    calculate_ecu_cycles()

<<<6000000>>>