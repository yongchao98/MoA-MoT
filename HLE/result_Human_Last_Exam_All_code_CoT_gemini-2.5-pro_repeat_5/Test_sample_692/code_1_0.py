def earliest_precolumbian_date():
    """
    This function provides information on the earliest known date recorded
    by a pre-Columbian civilization in the Americas.
    """
    # The earliest known calendar notation
    date_notation = "7 Deer"
    civilization = "Maya"
    location = "the Las Pinturas pyramid at San Bartolo, Guatemala"
    time_period = "3rd century BCE"

    # The earliest known Long Count date (for context)
    long_count_stela = "Stela 2 from Chiapa de Corzo, Mexico"
    long_count_date_mayan = "7.16.3.2.13"
    long_count_date_gregorian = "December 10, 36 BCE"

    print("The earliest known date recorded in an aboriginal American writing system is the '7 Deer' glyph.")
    print("\n--- Details of the Earliest Finding ---")
    print(f"Date Notation: {date_notation}")
    print(f"Civilization: {civilization}")
    print(f"Location: Found on painted mural fragments from {location}")
    print(f"Estimated Age: The finding is dated to the {time_period}, making it the oldest known calendar record from the Americas.")

    print("\n--- Previous Record Holder (for context) ---")
    print(f"For many years, one of the oldest known dates was a Long Count date from {long_count_stela}.")
    print(f"Mayan Long Count: {long_count_date_mayan[0]}.{long_count_date_mayan[2:5]}.{long_count_date_mayan[6]}.{long_count_date_mayan[8]}.{long_count_date_mayan[10:]}")
    print(f"Gregorian Equivalent: {long_count_date_gregorian}")

if __name__ == "__main__":
    earliest_precolumbian_date()