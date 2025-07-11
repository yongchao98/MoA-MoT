def find_town():
    """
    Identifies the correct town based on population and proximity constraints.

    The house in the image is located near Knox, Indiana. The prompt requires
    a town with a population over 5,000. Knox's population is too small.
    This function finds the closest nearby town that meets the population criteria.
    """
    # Data for nearby towns including their population (approx. 2020 census)
    # and approximate driving distance in miles from the house location.
    nearby_towns = {
        "Knox": {"population": 3650, "distance_miles": 4},
        "Plymouth": {"population": 10214, "distance_miles": 15},
        "La Porte": {"population": 22471, "distance_miles": 20},
        "Valparaiso": {"population": 34151, "distance_miles": 25}
    }

    # The minimum population required by the prompt.
    min_population = 5000

    # Filter the dictionary to find towns with a population greater than min_population.
    eligible_towns = {}
    for town, data in nearby_towns.items():
        if data["population"] > min_population:
            eligible_towns[town] = data

    # Find the closest town from the eligible list.
    # We are looking for the minimum distance among the eligible towns.
    if not eligible_towns:
        print("No nearby towns meet the population criteria.")
        return

    # Use a lambda function to find the key (town name) corresponding to the minimum distance.
    closest_town_name = min(eligible_towns, key=lambda town: eligible_towns[town]["distance_miles"])

    # This part fulfills the strange request to "output each number in the final equation"
    # by showing the numbers involved in the final selection.
    print(f"Checking eligible towns (population > {min_population}):")
    for town, data in eligible_towns.items():
        print(f"- {town}: population = {data['population']}, distance = {data['distance_miles']} miles")

    final_distance = eligible_towns[closest_town_name]['distance_miles']
    print(f"\nThe minimum distance is {final_distance} miles.")
    print("\nFinal Answer:")
    print(closest_town_name)

find_town()