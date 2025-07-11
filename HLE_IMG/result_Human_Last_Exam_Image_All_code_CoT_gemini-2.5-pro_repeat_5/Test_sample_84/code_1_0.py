def find_city_with_similar_coat_of_arms():
    """
    Identifies the city whose coat of arms is most similar to the one depicted.

    The image shows a coat of arms with the following key features on its shield:
    - Divided horizontally (per fess).
    - The top part is blue (azure) and contains a yellow/gold crescent moon.
    - The bottom part is yellow/gold (or) and contains a black crab.

    This specific design directly matches the coat of arms of a town in Croatia.
    """
    city_name = "Hrvatska Kostajnica"
    country = "Croatia"
    
    print(f"The coat of arms in the image is very similar to that of the city: {city_name}, {country}.")

find_city_with_similar_coat_of_arms()