def solve():
    """
    This function identifies the cities in the image and formats the output.
    """
    # The city where the picture was likely taken from, based on the vantage point.
    # The view is from the northwest of the background city, across the water.
    city_picture_taken_in = "Sausalito"

    # The city in the background, identifiable by the prominent illuminated tower
    # which strongly resembles Coit Tower.
    city_in_background = "San Francisco"

    # Print the result in the specified format: (city the picture was taken in, city where the background is)
    print(f"({city_picture_taken_in}, {city_in_background})")

solve()