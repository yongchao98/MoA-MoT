def calculate_ecu_cycles():
    """
    Calculates the maximum theoretical ECU cycles available between engine revolutions.
    """
    # Given values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450   # MegaHertz

    # Step 1: Convert Engine RPM to Revolutions per Second (RPS)
    # This also represents the number of interrupts per second.
    revolutions_per_second = engine_rpm / 60

    # Step 2: Convert ECU Speed from MHz to Hz (cycles per second)
    ecu_speed_hz = ecu_speed_mhz * 1_000_000

    # Step 3: Calculate the maximum cycles available per revolution (between interrupts)
    # Cycles per Revolution = (Total Cycles per Second) / (Revolutions per Second)
    max_cycles = ecu_speed_hz / revolutions_per_second

    # Print the explanation and the final equation with numbers
    print(f"The engine runs at {engine_rpm} RPM, and the ECU runs at {ecu_speed_mhz} MHz.")
    print("An interrupt happens once per revolution.")
    print("\nStep 1: Calculate revolutions per second (RPS)")
    print(f"{engine_rpm} RPM / 60 seconds/minute = {revolutions_per_second} RPS")
    
    print("\nStep 2: Calculate ECU cycles per second (Hz)")
    print(f"{ecu_speed_mhz} MHz * 1,000,000 Hz/MHz = {int(ecu_speed_hz)} Hz")

    print("\nStep 3: Calculate the maximum cycles available per revolution")
    print("The formula is: Cycles per Revolution = ECU Speed (Hz) / Engine Speed (RPS)")
    
    # Final equation with numbers
    print("\nFinal Equation:")
    print(f"Maximum Cycles = {int(ecu_speed_hz)} / {int(revolutions_per_second)}")

    print(f"\nResult:")
    print(f"The theoretical maximum number of cycles is {int(max_cycles)}")

# Execute the function
calculate_ecu_cycles()
<<<6000000>>>