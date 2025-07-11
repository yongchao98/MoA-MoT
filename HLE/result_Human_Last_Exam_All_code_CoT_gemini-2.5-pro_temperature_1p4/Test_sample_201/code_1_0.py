import pandas as pd

def solve():
    """
    Analyzes pond characteristics to determine the greatest mosquito threat.
    """
    # Create a list of dictionaries to represent the ponds
    ponds_data = [
        {'id': 'A', 'size_ft': 10, 'age_yrs': 1},
        {'id': 'C', 'size_ft': 30, 'age_yrs': 1},
        {'id': 'D', 'size_ft': 10, 'age_yrs': 5},
        {'id': 'E', 'size_ft': 30, 'age_yrs': 5}
    ]

    # Initialize variables to find the pond with the maximum threat
    max_threat_score = -1
    best_pond = None
    winning_equation = ""

    print("Analyzing ponds to find the greatest medical threat...")
    print("Threat is proportional to both surface area and age of the pond.\n")

    # Calculate a threat score for each pond
    for pond in ponds_data:
        area = pond['size_ft'] * pond['size_ft']
        age = pond['age_yrs']
        # Threat score is a product of area and age, representing combined favorability
        threat_score = area * age
        
        # Check if this pond has the highest score so far
        if threat_score > max_threat_score:
            max_threat_score = threat_score
            best_pond = pond
            winning_equation = f"Threat Score = Area * Age = {pond['size_ft']} ft * {pond['size_ft']} ft * {pond['age_yrs']} years = {threat_score}"


    # Print the result
    print("Conclusion:")
    print(f"The pond representing the greatest medical threat is Pond {best_pond['id']}.")
    print(f"It is {best_pond['size_ft']} feet square and {best_pond['age_yrs']} years old.")
    print("This pond has the largest surface area and is the oldest, allowing for the most established insect community.")
    print("\nCalculation for the winning pond:")
    print(winning_equation)


solve()
<<<E>>>