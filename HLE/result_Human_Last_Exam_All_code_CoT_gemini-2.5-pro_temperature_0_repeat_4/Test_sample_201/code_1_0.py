def calculate_threat():
    """
    Calculates a relative threat score for different ponds based on size and age.
    """

    # Define the options as a dictionary
    # Format: { 'Option': (side_length, age_in_years) }
    options = {
        'A': (10, 1),
        'C': (30, 1),
        'D': (10, 5),
        'E': (30, 5)
    }

    # Assign a multiplier for age. Older ponds are more established.
    # Let's use a simple factor where a 5-year-old pond is twice as established as a 1-year-old.
    age_factors = {1: 1, 5: 2}

    print("Calculating a 'Threat Score' for each pond.")
    print("Threat Score = (Side * Side) * Age_Factor\n")

    threat_scores = {}
    for option, (side, age) in options.items():
        surface_area = side * side
        age_factor = age_factors[age]
        score = surface_area * age_factor
        threat_scores[option] = score
        print(f"Pond {option} ({side}ft square, {age} year(s) old):")
        print(f"Threat Score = ({side} * {side}) * {age_factor} = {score}\n")

    # Find the option with the highest score
    max_threat_option = max(threat_scores, key=threat_scores.get)
    
    print(f"Conclusion: Pond {max_threat_option} has the highest threat score, indicating it would host the highest abundance of mosquitoes.")

calculate_threat()
<<<E>>>