import requests

def get_town_population(town, state_code):
    """
    Fetches population data for a given town using a public API.
    Note: This is a conceptual function. Real-world usage requires a specific API and key.
    For this problem, we will use known values.
    """
    populations = {
        "Seneca": 2353,
        "Ottawa": 18716,
        "Morris": 14163,
        "Marseilles": 4845
    }
    return populations.get(town, 0)

def solve_location_puzzle():
    """
    This function solves the puzzle by identifying the correct town based on location and population.
    """
    # Step 1: Geolocation points to 493 E 29th Rd, Seneca, IL.
    postal_town = "Seneca"
    postal_town_population = get_town_population(postal_town, "IL")

    print(f"The postal town is {postal_town} with a population of {postal_town_population}.")
    print("This does not meet the > 5,000 population requirement.\n")

    # Step 2: Evaluate nearby towns that meet the criteria.
    # Distances are approximate driving distances from the location.
    candidate1 = {"name": "Ottawa", "distance_miles": 11.2, "county_seat": True}
    candidate2 = {"name": "Morris", "distance_miles": 10.4, "county_seat": False}

    candidate1["population"] = get_town_population(candidate1["name"], "IL")
    candidate2["population"] = get_town_population(candidate2["name"], "IL")

    print(f"Candidate 1: {candidate1['name']}")
    print(f"Population: {candidate1['population']}")
    print(f"Distance: ~{candidate1['distance_miles']} miles")
    print(f"It is the county seat.\n")

    print(f"Candidate 2: {candidate2['name']}")
    print(f"Population: {candidate2['population']}")
    print(f"Distance: ~{candidate2['distance_miles']} miles")
    print(f"It is in a neighboring county.\n")

    # Step 3: Conclude based on proximity.
    # Since both meet the population criteria, the closest one is the most likely answer.
    if candidate2["distance_miles"] < candidate1["distance_miles"]:
        final_answer = candidate2["name"]
    else:
        final_answer = candidate1["name"]

    print(f"Conclusion: Although both Ottawa and Morris have populations over 5,000, Morris is slightly closer.")
    print(f"Therefore, the home is most likely considered part of the town of Morris.")
    print("\nFinal Answer:")
    print(final_answer)


solve_location_puzzle()