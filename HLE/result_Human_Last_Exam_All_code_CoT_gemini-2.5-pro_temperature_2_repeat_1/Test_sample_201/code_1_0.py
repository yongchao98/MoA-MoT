import sys

def solve_mosquito_threat():
    """
    Calculates the 'threat score' for different ponds to identify the one
    posing the greatest medical threat from mosquitoes.
    """

    # Pond configurations from the answer choices
    # Each dictionary contains an ID, side length in feet, and age in years.
    ponds = [
        {'id': 'A', 'side': 10, 'age': 1},
        {'id': 'C', 'side': 30, 'age': 1},
        {'id': 'D', 'side': 10, 'age': 5},
        {'id': 'E', 'side': 30, 'age': 5},
    ]

    max_threat_score = -1
    best_pond = None

    print("Calculating the medical threat score for each pond...\n")

    for pond in ponds:
        # Surface area is side * side for a square pond.
        surface_area = pond['side'] * pond['side']
        
        # Age represents the maturity of the insect community. A higher age
        # implies a more established community, thus a higher multiplier.
        age_factor = pond['age']

        # The threat score is a product of the area available for larvae and the
        # maturity of the ecosystem.
        threat_score = surface_area * age_factor
        
        print(f"Pond {pond['id']} ({pond['side']}ft square, {pond['age']} year(s) old):")
        print(f"  Threat Score = Surface Area * Age Factor")
        print(f"  Threat Score = {surface_area} * {age_factor} = {threat_score}\n")

        # Keep track of the pond with the highest score
        if threat_score > max_threat_score:
            max_threat_score = threat_score
            best_pond = pond
    
    print("-" * 30)
    print(f"Conclusion: Pond {best_pond['id']} has the highest threat score.")
    print(f"This pond is {best_pond['side']} feet square and {best_pond['age']} years old, offering the largest surface area")
    print("and the most established insect community, leading to the greatest medical threat.")
    
    # The final answer is wrapped in '<<<' and '>>>'
    # sys.stdout.write needs to be used to avoid extra newlines from print()
    sys.stdout.write(f"<<<{best_pond['id']}>>>")

solve_mosquito_threat()