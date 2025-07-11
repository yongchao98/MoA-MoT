import sys

def solve_mosquito_threat():
    """
    This script determines which pond represents the greatest medical threat
    by identifying the one likely to host the highest mosquito abundance.

    Methodology:
    1. Mosquito abundance is primarily influenced by available habitat (surface area)
       and the maturity of the aquatic ecosystem (age).
    2. A larger surface area (side * side) can support more mosquito larvae.
    3. An older pond has a more established ecosystem with more food for larvae.
    4. A 'threat score' is calculated for each pond by multiplying its surface area by its age.
       Threat Score = (side * side) * age
    5. The pond with the highest score is the one with the highest threat.
    """

    ponds = {
        'A': {'side': 10, 'age': 1, 'description': "10 feet square, one year old"},
        'C': {'side': 30, 'age': 1, 'description': "30 feet square, one year old"},
        'D': {'side': 10, 'age': 5, 'description': "10 feet square, five years old"},
        'E': {'side': 30, 'age': 5, 'description': "30 feet square, five years old"}
    }

    max_threat_score = -1
    best_pond_id = None

    print("Calculating threat score for each pond (Area * Age):")

    for pond_id, properties in ponds.items():
        side = properties['side']
        age = properties['age']
        area = side * side
        threat_score = area * age
        print(f"- Pond {pond_id} ({properties['description']}): Score = {area} * {age} = {threat_score}")

        if threat_score > max_threat_score:
            max_threat_score = threat_score
            best_pond_id = pond_id

    print("\n--- Conclusion ---")
    best_pond = ponds[best_pond_id]
    side = best_pond['side']
    age = best_pond['age']

    print(f"The pond representing the greatest medical threat is Pond {best_pond_id}.")
    print("This is because it has the largest surface area and is the oldest, allowing for the most established insect community.")
    print(f"The winning calculation is for the {best_pond['description']} pond:")
    # The final print statement showing each number in the final equation
    print(f"Threat Score = {side} feet * {side} feet * {age} years = {max_threat_score}")

solve_mosquito_threat()
sys.stdout.write("<<<E>>>")