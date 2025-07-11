import sys

def solve():
    """
    This script calculates a "Threat Score" for each pond to determine which
    poses the greatest risk. The score is based on the principle that larger
    and older ponds are better habitats for mosquitoes.
    """

    # Pond data from the answer choices
    ponds = {
        'A': {'size': 10, 'age': 1},
        'C': {'size': 30, 'age': 1},
        'D': {'size': 10, 'age': 5},
        'E': {'size': 30, 'age': 5},
    }

    # A simple scoring system.
    # For size, we use a relative score: 1 for 10ft, 3 for 30ft.
    # For age, we use the age in years as the score.
    # The total score is the product, as the effects are cumulative.
    print("Calculating a Threat Score for each pond.")
    print("The score is based on the equation: Size_Factor * Age_Factor = Total_Score\n")

    scores = {}
    best_pond = None
    max_score = -1

    for key, properties in ponds.items():
        size = properties['size']
        age = properties['age']

        # Assigning factors
        size_factor = 3 if size == 30 else 1
        age_factor = age

        # The final equation and result
        total_score = size_factor * age_factor
        scores[key] = total_score
        
        print(f"Pond {key} ({size} feet square, {age} year(s) old):")
        # As requested, printing each number in the final equation
        print(f"  Equation: {size_factor} * {age_factor} = {total_score}")

        if total_score > max_score:
            max_score = total_score
            best_pond = key
            
    print(f"\nThe pond with the highest score is '{best_pond}'.")
    print("This pond (30 feet square, five years old) represents the greatest medical threat by providing the largest and most established habitat for mosquitoes.")

solve()