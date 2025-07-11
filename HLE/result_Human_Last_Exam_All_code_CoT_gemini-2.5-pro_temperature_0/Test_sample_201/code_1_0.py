def calculate_threat():
    """
    Calculates a threat score for different ponds based on size and age
    to determine which poses the greatest medical threat from mosquitoes.
    """

    # Define the characteristics of the ponds from the answer choices
    ponds = {
        'A': {'side_length': 10, 'age': 1},
        'C': {'side_length': 30, 'age': 1},
        'D': {'side_length': 10, 'age': 5},
        'E': {'side_length': 30, 'age': 5},
    }

    print("Calculating a 'threat score' for each pond.")
    print("Threat Score = Surface Area * Age\n")

    max_threat_score = 0
    highest_threat_pond = None

    # Calculate and print the threat score for each pond
    for key, properties in ponds.items():
        side = properties['side_length']
        age = properties['age']
        
        # Surface area of a square pond is side * side
        surface_area = side * side
        
        # Threat score is a product of area and age
        threat_score = surface_area * age
        
        print(f"Pond {key}:")
        print(f"  Equation: {surface_area} (Area) * {age} (Age) = {threat_score} (Threat Score)")
        
        if threat_score > max_threat_score:
            max_threat_score = threat_score
            highest_threat_pond = key

    print("\nConclusion:")
    print(f"Pond {highest_threat_pond} has the highest threat score of {max_threat_score}.")
    print("This is because it is both the largest and the oldest, providing the most suitable environment for a large, established mosquito population.")

calculate_threat()