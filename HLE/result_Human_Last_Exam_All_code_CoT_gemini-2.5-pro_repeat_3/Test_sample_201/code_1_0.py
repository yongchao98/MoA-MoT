import operator

def calculate_threat():
    """
    Calculates a threat score for different ponds based on size and age.
    """
    ponds = {
        'A': {'side': 10, 'age': 1, 'description': '10 feet square, one year old'},
        'C': {'side': 30, 'age': 1, 'description': '30 feet square, one year old'},
        'D': {'side': 10, 'age': 5, 'description': '10 feet square, five years old'},
        'E': {'side': 30, 'age': 5, 'description': '30 feet square, five years old'},
    }

    threat_scores = {}

    print("Calculating the relative medical threat for each pond.")
    print("Threat is modeled as: (Side * Side) * Age_Factor\n")

    for key, properties in ponds.items():
        side = properties['side']
        age = properties['age']
        
        # Calculate surface area
        area = side * side
        # Calculate threat score
        threat_score = area * age
        
        threat_scores[key] = threat_score
        
        print(f"Pond {key} ({properties['description']}):")
        # The final equation with each number
        print(f"Threat Score = {side} * {side} * {age} = {threat_score}\n")

    # Find the pond with the highest threat score
    max_threat_pond = max(threat_scores.items(), key=operator.itemgetter(1))
    
    print(f"The pond with the highest threat score is Pond {max_threat_pond[0]}.")
    print("This is because it has the largest surface area for breeding and is the oldest, allowing for a more established insect community.")

calculate_threat()