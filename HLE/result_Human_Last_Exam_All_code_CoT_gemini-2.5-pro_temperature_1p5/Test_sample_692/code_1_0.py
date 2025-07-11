import textwrap

def find_earliest_date():
    """
    This function presents information about the earliest known date recorded
    by a pre-Columbian civilization in an aboriginal writing system.
    """
    civilization = "Maya"
    site = "San Bartolo, Guatemala"
    date_notation = "7 Deer"
    calendar_system = "260-day Tzolk'in (sacred calendar)"
    approximate_date = "c. 250 BC"
    context = (
        "This discovery, announced in April 2022, was made on painted plaster "
        "fragments from the Las Pinturas pyramid at the San Bartolo archaeological site. "
        "Radiocarbon dating of the charcoal associated with the fragments places them "
        "to between 300 and 200 BC, making the '7 Deer' notation the earliest "
        "secure evidence of a calendar being recorded in the Maya region."
    )

    # For context, this discovery predates the earliest known Long Count date.
    previous_record_holder = "Chiapa de Corzo Stela 2"
    previous_record_date_long_count = "7.16.3.2.13"
    previous_record_date_gregorian = "December 10, 36 BC"

    print("--- Earliest Known Pre-Columbian Recorded Date ---")
    print(f"Civilization: {civilization}")
    print(f"Archaeological Site: {site}")
    print(f"Recorded Date Notation: '{date_notation}'")
    print(f"Calendar System: {calendar_system}")
    print(f"Approximate Gregorian Date: {approximate_date}")

    print("\nContext:")
    print(textwrap.fill(context, 70))

    print("\n--- Previous Record Holder (Earliest Long Count Date) ---")
    print(f"Artifact: {previous_record_holder}")
    print(f"Mayan Long Count Date: {previous_record_date_long_count}")
    print(f"Equivalent Gregorian Date: {previous_record_date_gregorian}")


if __name__ == "__main__":
    find_earliest_date()
