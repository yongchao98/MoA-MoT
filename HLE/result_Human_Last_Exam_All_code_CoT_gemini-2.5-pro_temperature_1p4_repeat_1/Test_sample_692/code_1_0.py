import sys

def find_earliest_date():
    """
    This function identifies the earliest known recorded date from a pre-Columbian
    civilization by comparing several key archaeological finds.
    """
    # Data on key archaeological finds with early recorded dates.
    # BCE dates are stored as negative numbers for easy chronological sorting.
    # The higher the absolute value, the earlier the date.
    artifacts = [
        {
            "site": "Tres Zapotes, Mexico",
            "artifact": "Stela C",
            "culture": "Epi-Olmec",
            "year_bce": 32,
            "date_notation": "Long Count date 7.16.6.16.18",
            "note": "For a long time, this was considered one of the oldest Long Count dates."
        },
        {
            "site": "Chiapa de Corzo, Mexico",
            "artifact": "Stela 2",
            "culture": "Epi-Olmec",
            "year_bce": 36,
            "date_notation": "Long Count date 7.16.3.2.13",
            "note": "Holds the record for the oldest known date from the Long Count calendar system."
        },
        {
            "site": "San Bartolo, Guatemala",
            "artifact": "Mural Fragments",
            "culture": "Maya",
            "year_bce": 250, # This is an approximate date from radiocarbon dating
            "date_notation": "Day '7 Deer' from the 260-day ritual calendar",
            "note": "Represents the earliest definitive evidence of a Mesoamerican calendar date, published in 2022."
        }
    ]

    # Find the artifact with the largest BCE year (i.e., the earliest date)
    earliest_artifact = max(artifacts, key=lambda x: x['year_bce'])

    site = earliest_artifact["site"]
    date_notation = earliest_artifact["date_notation"]
    year = earliest_artifact["year_bce"]

    # Extracting numbers from the date notation for the final output string.
    # This interprets the user's request to "output each number in the final equation".
    numbers_in_date = [s for s in date_notation.split() if s.isdigit()]

    print(f"The earliest known date recorded by a pre-Columbian civilization was found at {site}.")
    print(f"The inscription is a date from the 260-day ritual calendar.")
    print("\n--- Date Details ---")
    print(f"Recorded Date: '{date_notation}'")
    print(f"Estimated Year: Approximately {year} BCE")
    # Fulfills the requirement to output the numbers in the "equation" or notation.
    print(f"The key number from the notation is: {numbers_in_date[0]}")
    print(f"The year associated with the finding is: {year}")


find_earliest_date()