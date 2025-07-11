def solve_city_puzzle():
    """
    Identifies the city where the picture was taken and the city in the background.
    """
    # Based on the landmark Coit Tower, the background city is San Francisco.
    city_in_background = "San Francisco"

    # The vantage point across the bay, consistent with the view, is Sausalito.
    city_taken_in = "Sausalito"

    # Print the result in the specified format.
    print(f"({city_taken_in}, {city_in_background})")

solve_city_puzzle()