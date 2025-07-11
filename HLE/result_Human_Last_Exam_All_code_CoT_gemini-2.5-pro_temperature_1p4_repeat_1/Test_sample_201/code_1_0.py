import sys

def solve():
    """
    This script evaluates which pond represents the greatest medical threat
    by identifying the one most likely to host the highest mosquito abundance.
    """

    # Define the ponds based on the answer choices
    ponds = {
        'A': {'size': 10, 'age': 1},
        'C': {'size': 30, 'age': 1},
        'D': {'size': 10, 'age': 5},
        'E': {'size': 30, 'age': 5}
    }

    # Initialize variables to find the best pond
    best_pond_label = ''
    max_threat_score = -1

    print("Analyzing factors contributing to mosquito abundance:")
    print("1. Habitat Size: A larger pond (30ft) provides more breeding ground than a smaller one (10ft).")
    print("2. Ecosystem Maturity: An older pond (5 years) has a more established insect community than a younger one (1 year).")
    print("\nTo find the greatest threat, we look for the pond that is both the largest AND the oldest.")
    print("We can model this by calculating a 'Threat Score' = size * age.\n")

    # Calculate threat score for each pond
    for label, properties in ponds.items():
        size = properties['size']
        age = properties['age']
        threat_score = size * age
        print(f"Pond {label} ({size}ft, {age}yr): Threat Score = {size} * {age} = {threat_score}")

        if threat_score > max_threat_score:
            max_threat_score = threat_score
            best_pond_label = label

    best_pond = ponds[best_pond_label]

    # Print the conclusion
    print(f"\nThe pond with the highest score is '{best_pond_label}'.")
    print(f"\nFinal Calculation for the Highest Threat Pond:")
    # The final equation output as requested
    final_equation_str = f"{best_pond['size']} feet square * {best_pond['age']} years old = {max_threat_score} Threat Score"
    print(final_equation_str)
    print("\nThis pond represents the greatest medical threat due to having the largest size and the most established ecosystem.")


solve()
<<<E>>>