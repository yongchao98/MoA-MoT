def find_artist_by_date():
    """
    This function identifies and prints the musical artist who performed at the
    Eccles Theater on a specific date.
    """
    venue_name = "Eccles Theater"
    month = "October"
    day = 12
    year = 2023
    artist_name = "The National"

    print(f"The musical artist that played at the {venue_name} on {month} {day}, {year} was: {artist_name}")

if __name__ == "__main__":
    find_artist_by_date()