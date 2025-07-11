import operator

def solve_mosquito_problem():
    """
    Calculates a threat score for different ponds to determine which
    poses the greatest medical threat from mosquitoes.
    """
    ponds = {
        'A': {'size_ft': 10, 'age_years': 1, 'description': '10 feet square, one year old.'},
        'C': {'size_ft': 30, 'age_years': 1, 'description': '30 feet square, one year old.'},
        'D': {'size_ft': 10, 'age_years': 5, 'description': '10 feet square, five years old.'},
        'E': {'size_ft': 30, 'age_years': 5, 'description': '30 feet square, five years old.'}
    }

    threat_scores = {}

    print("Calculating threat scores for each pond...")
    print("Formula: Threat Score = (Side * Side) * Age\n")

    for key, properties in ponds.items():
        side = properties['size_ft']
        age = properties['age_years']
        area = side * side
        score = area * age
        threat_scores[key] = score
        print(f"Pond {key} ({properties['description']}):")
        print(f"Threat Score = ({side} * {side}) * {age} = {score}\n")

    # Find the pond with the highest score
    best_option_key = max(threat_scores.items(), key=operator.itemgetter(1))[0]
    best_option_description = ponds[best_option_key]['description']

    print(f"Conclusion: The pond with the highest threat score is Option {best_option_key}.")
    print(f"This pond is {best_option_description}, as it is both the largest and the oldest, providing the best environment for a large mosquito population.")

solve_mosquito_problem()
<<<E>>>