def analyze_mosquito_threat():
    """
    Analyzes pond options to determine which poses the greatest medical threat
    based on mosquito abundance.
    """
    # Assign scores based on the factors. Larger and older are better for mosquitoes.
    # A 30ft pond is 9x larger in area (900 sq ft vs 100 sq ft), so we'll give it a higher score.
    # An older pond is explicitly stated to be more established.
    size_scores = {10: 1, 30: 9}
    age_scores = {1: 1, 5: 5}

    # Define the pond options from the question
    options = {
        'A': {'size': 10, 'age': 1},
        'C': {'size': 30, 'age': 1},
        'D': {'size': 10, 'age': 5},
        'E': {'size': 30, 'age': 5}
    }

    results = {}
    print("Calculating a 'Threat Score' for each pond.")
    print("Threat Score = (Size Score) * (Age Score)\n")

    # Calculate the threat score for each option
    for key, props in options.items():
        size = props['size']
        age = props['age']
        size_score = size_scores[size]
        age_score = age_scores[age]
        total_score = size_score * age_score
        
        # Store the result and the equation string
        results[key] = total_score
        print(f"Option {key} ({size}ft, {age}yr):")
        # The final code needs to output each number in the final equation
        print(f"  Equation: {size_score} * {age_score} = {total_score}")

    # Find the option with the highest score
    best_option = max(results, key=results.get)

    print(f"\nConclusion: Option {best_option} has the highest threat score.")
    print("A pond that is both larger (30 feet square) and older (five years old) will have the most established insect community and the largest habitat, supporting the highest abundance of mosquitoes.")

analyze_mosquito_threat()