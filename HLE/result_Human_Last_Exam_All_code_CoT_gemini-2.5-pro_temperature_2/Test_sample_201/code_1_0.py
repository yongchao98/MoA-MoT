import math

def analyze_mosquito_threat():
    """
    This script analyzes different pond characteristics to determine which
    poses the greatest medical threat based on potential mosquito abundance.

    The model assumes that threat is proportional to abundance, which is influenced by:
    1. Surface Area: Larger ponds support more larvae.
    2. Age: Older ponds have more established insect communities, making them
       better habitats.
    """

    ponds = {
        'A': {'size': 10, 'age': 1, 'description': "10 feet square, one year old."},
        'C': {'size': 30, 'age': 1, 'description': "30 feet square, one year old."},
        'D': {'size': 10, 'age': 5, 'description': "10 feet square, five years old."},
        'E': {'size': 30, 'age': 5, 'description': "30 feet square, five years old."}
    }

    # To model the "more established communities", we'll use a simple age factor.
    # Let's assume a 5-year-old pond is twice as productive as a 1-year-old pond.
    age_factor_map = {1: 1, 5: 2}

    highest_threat_score = -1
    best_option = None

    print("Calculating Threat Score for each pond...")
    print("Formula: Threat Score = (Surface Area) * (Age Factor)")
    print("-" * 50)

    for option, details in ponds.items():
        size = details['size']
        age = details['age']
        description = details['description']

        surface_area = size ** 2
        age_factor = age_factor_map[age]
        threat_score = surface_area * age_factor

        print(f"Pond Option {option}: {description}")
        # Here we output the numbers in the final equation as requested.
        print(f"  Equation: {surface_area} (from {size}x{size}) * {age_factor} (for age {age}) = {threat_score}")

        if threat_score > highest_threat_score:
            highest_threat_score = threat_score
            best_option = option

    print("-" * 50)
    print(f"The pond representing the greatest medical threat is Option {best_option},")
    print(f"as it has the largest size and is older, resulting in the highest threat score.")

analyze_mosquito_threat()