import math

def calculate_threat_score():
    """
    Calculates a threat score for different ponds to determine which
    poses the greatest risk for mosquito abundance.
    """
    # Define the characteristics of the ponds from the answer choices
    ponds = [
        {"option": "A", "side_length_ft": 10, "age_years": 1},
        {"option": "C", "side_length_ft": 30, "age_years": 1},
        {"option": "D", "side_length_ft": 10, "age_years": 5},
        {"option": "E", "side_length_ft": 30, "age_years": 5},
    ]

    highest_score = -1
    best_pond = None

    print("Calculating threat scores for each pond...\n")

    for pond in ponds:
        # Calculate surface area (side_length * side_length)
        area = pond["side_length_ft"] * pond["side_length_ft"]
        age = pond["age_years"]

        # The threat is proportional to both area and age.
        # We can model this with a simple score: Threat = Area * Age
        threat_score = area * age
        
        # Keep track of the pond with the highest score
        if threat_score > highest_score:
            highest_score = threat_score
            best_pond = pond

        print(f"Pond {pond['option']}:")
        print(f"  - Characteristics: {pond['side_length_ft']} feet square, {pond['age_years']} year(s) old.")
        print(f"  - Threat Score Calculation: {area} (Area) * {age} (Age) = {threat_score}")
        print("-" * 20)

    print(f"\nConclusion:")
    print(f"The pond with the highest threat score is Option {best_pond['option']}.")
    print(f"This is the {best_pond['side_length_ft']} feet square, {best_pond['age_years']} years old pond, as it combines the largest area with the oldest, most established ecosystem, supporting the highest mosquito abundance.")

if __name__ == '__main__':
    calculate_threat_score()
<<<E>>>