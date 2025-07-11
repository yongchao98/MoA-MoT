def calculate_rarest_noble_gas_percentage():
    """
    This script identifies the rarest stable noble gas and calculates its
    percentage in Earth's atmosphere.
    """
    gas_name = "Xenon"
    abundance_ppm = 0.087  # Abundance in parts per million

    # Convert ppm to percentage
    # percentage = (parts per million / 1,000,000) * 100
    divisor = 1000000
    multiplier = 100
    percentage = (abundance_ppm / divisor) * multiplier

    print(f"The rarest stable noble gas on Earth is {gas_name}.")
    print(f"Its approximate abundance in the atmosphere is {abundance_ppm} parts per million (ppm).")
    print("This value is used as a standard measure for its presence in terrestrial matter.")
    print("\nThe calculation to convert its abundance from ppm to a percentage is shown below:")
    # The format below prints each number in the final equation, as requested.
    print(f"{abundance_ppm} / {divisor} * {multiplier} = {percentage}%")

calculate_rarest_noble_gas_percentage()