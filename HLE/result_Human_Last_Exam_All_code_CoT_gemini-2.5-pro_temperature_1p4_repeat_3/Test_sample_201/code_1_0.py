def calculate_mosquito_threat():
    """
    Calculates a threat score for different ponds based on size and age
    to determine which poses the greatest medical threat.
    """

    # Pond data based on the answer choices
    # Format: {ID: {"size_ft": value, "age_years": value}}
    ponds = {
        'A': {'size_ft': 10, 'age_years': 1},
        'C': {'size_ft': 30, 'age_years': 1},
        'D': {'size_ft': 10, 'age_years': 5},
        'E': {'size_ft': 30, 'age_years': 5},
    }

    best_pond = None
    max_threat_score = -1
    final_equation_str = ""

    print("Analyzing mosquito threat based on pond characteristics.")
    print("Threat Score = Area (sq ft) * Age Factor (years)")
    print("--------------------------------------------------")

    # Iterate through each pond to calculate its threat score
    for key, data in ponds.items():
        size = data['size_ft']
        age = data['age_years']
        
        # Area is a primary factor for habitat size
        area = size * size
        
        # Age is a multiplier, as older ponds have more established communities
        # The age itself can serve as the factor
        threat_score = area * age
        
        # Print the calculation for each pond
        equation_str = f"Pond {key} ({size}ft, {age}yr): Threat = {area} * {age} = {threat_score}"
        print(equation_str)

        if threat_score > max_threat_score:
            max_threat_score = threat_score
            best_pond = key
            final_equation_str = equation_str

    print("--------------------------------------------------")
    print(f"Conclusion: Pond {best_pond} has the highest threat score.")
    print("This is because it is both the largest and the oldest, providing the most ideal conditions for a large mosquito population.")
    print("\nFinal equation for the highest threat pond:")
    # The requirement is to output the numbers in the final equation.
    # The string `final_equation_str` already contains this formatted information.
    print(final_equation_str)


calculate_mosquito_threat()
<<<E>>>