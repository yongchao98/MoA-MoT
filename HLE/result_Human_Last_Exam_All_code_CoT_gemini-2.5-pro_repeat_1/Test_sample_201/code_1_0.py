import operator

def solve_mosquito_threat():
    """
    Calculates the medical threat level for different ponds based on size and age.
    """
    # Define the ponds based on the answer choices.
    # Format: {'Choice': ('Description', {'size': value, 'age': value})}
    ponds = {
        'A': ('10 feet square, one year old', {'size': 10, 'age': 1}),
        'C': ('30 feet square, one year old', {'size': 30, 'age': 1}),
        'D': ('10 feet square, five years old', {'size': 10, 'age': 5}),
        'E': ('30 feet square, five years old', {'size': 30, 'age': 5})
    }

    threat_scores = {}

    print("Calculating threat scores for each pond...")
    print("Threat Score = Pond Size (ft) * Pond Age (years)\n")

    # Calculate the threat score for each pond
    for choice, (description, attributes) in ponds.items():
        size = attributes['size']
        age = attributes['age']
        score = size * age
        threat_scores[choice] = score
        # The final code needs to output each number in the final equation
        print(f"Pond {choice} ({description}):")
        print(f"  Calculation: {size} * {age} = {score}")
        print("-" * 20)

    # Find the pond with the highest threat score
    if not threat_scores:
        print("No ponds to evaluate.")
        return

    # Find the key (pond choice) with the maximum value (score)
    max_threat_pond = max(threat_scores.items(), key=operator.itemgetter(1))[0]

    print(f"\nConclusion:")
    print(f"The pond with the highest threat score is Pond {max_threat_pond}.")
    print("A larger and older pond provides a more extensive and stable habitat,")
    print("leading to a higher abundance of mosquitoes and thus a greater medical threat.")


solve_mosquito_threat()