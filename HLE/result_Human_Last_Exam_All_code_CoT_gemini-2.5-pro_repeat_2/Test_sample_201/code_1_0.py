def calculate_threat():
    """
    Calculates and compares the medical threat level of different ponds
    based on their size and age.
    """
    # Pond characteristics based on answer choices
    # Format: [choice_letter, side_length, age_in_years, description]
    ponds = [
        ['A', 10, 1, '10 feet square, one year old'],
        ['C', 30, 1, '30 feet square, one year old'],
        ['D', 10, 5, '10 feet square, five years old'],
        ['E', 30, 5, '30 feet square, five years old']
    ]

    max_threat = -1
    best_pond_choice = ''
    
    print("Calculating a 'threat score' for each pond based on the formula: Area * Age.\n")

    for choice, side, age, description in ponds:
        area = side * side
        threat_score = area * age
        
        print(f"Pond {choice} ({description}):")
        # The request is to output each number in the final equation
        print(f"Threat Score = Area * Age = {side} * {side} * {age} = {threat_score}\n")

        if threat_score > max_threat:
            max_threat = threat_score
            best_pond_choice = choice

    print("Conclusion:")
    print(f"The pond with the highest threat score is Pond {best_pond_choice}, which is the 30 feet square, five-year-old pond.")
    print("This is because a larger area and an older, more established ecosystem support a higher abundance of mosquitoes.")

calculate_threat()