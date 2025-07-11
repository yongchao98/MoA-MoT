def calculate_ecu_cycles():
    """
    Calculates the theoretical maximum number of CPU cycles per engine revolution.
    """
    # 1. Establish Given Values
    engine_rpm = 4500  # Revolutions Per Minute
    ecu_speed_mhz = 450  # Megahertz

    # 2. Standardize Units
    # Convert engine speed to Revolutions Per Second (RPS)
    seconds_in_minute = 60
    engine_rps = engine_rpm / seconds_in_minute

    # Convert ECU speed to cycles per second (Hertz)
    hertz_in_megahertz = 1_000_000
    ecu_speed_hz = ecu_speed_mhz * hertz_in_megahertz

    # 3. Calculate Cycles per Revolution
    # This is the number of CPU cycles divided by the number of revolutions in the same time period (one second).
    max_cycles_per_revolution = ecu_speed_hz / engine_rps

    # 4. Display the Result
    # Print the equation with all the numbers used in the calculation.
    print(f"Engine Speed: {engine_rpm} RPM")
    print(f"ECU Processor Speed: {ecu_speed_mhz} MHz")
    print("\nCalculation:")
    print(f"({int(ecu_speed_hz)} cycles/sec) / ({int(engine_rpm)} rev/min / {seconds_in_minute} sec/min) = {int(max_cycles_per_revolution)} cycles/revolution")

    # The final numerical answer is part of the string above,
    # but we can print it separately if needed.
    # print(f"\nFinal Answer: {int(max_cycles_per_revolution)}")


if __name__ == '__main__':
    calculate_ecu_cycles()
    # The result of the calculation is 6,000,000
    print("\n<<<6000000>>>")
