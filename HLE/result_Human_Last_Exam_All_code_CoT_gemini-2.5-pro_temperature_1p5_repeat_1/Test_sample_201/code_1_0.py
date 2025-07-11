import operator

def calculate_threat():
    """
    Calculates a "Threat Score" for different ponds to determine which
    poses the greatest medical risk based on mosquito abundance.
    
    The threat is modeled as a product of surface area and an age factor.
    """
    
    # Pond data: [description, side_length (ft), age (years)]
    ponds = {
        "A": ["10 feet square, one year old", 10, 1],
        "C": ["30 feet square, one year old", 30, 1],
        "D": ["10 feet square, five years old", 10, 5],
        "E": ["30 feet square, five years old", 30, 5],
    }
    
    # We model the "established community" aspect with an age factor.
    # A 5-year-old pond is more established than a 1-year-old one.
    # Let's assign a simple multiplier: Age 1 -> factor 1, Age 5 -> factor 2
    age_factors = {1: 1, 5: 2}
    
    results = {}

    print("Calculating Threat Score for each pond...\nThreat Score = (Side * Side) * Age Factor\n")

    for key, data in ponds.items():
        description, side, age = data
        area = side * side
        age_factor = age_factors[age]
        
        # Calculate the threat score
        threat_score = area * age_factor
        results[key] = threat_score
        
        print(f"Pond {key} ({description}):")
        print(f"Calculation: ({side} * {side}) * {age_factor} = {threat_score}")
        print("-" * 20)

    # Find the pond with the highest score
    highest_threat_pond = max(results.items(), key=operator.itemgetter(1))
    
    print(f"\nConclusion: Pond {highest_threat_pond[0]} has the highest threat score ({highest_threat_pond[1]}).")
    print("This is because it has both the largest surface area and is older, allowing for a more established insect community.")

calculate_threat()