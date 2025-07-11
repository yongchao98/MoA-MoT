def solve_mosquito_threat():
    """
    Calculates a "threat score" for different ponds to determine which
    poses the greatest medical threat based on mosquito abundance.
    The score is based on the principle that larger and older ponds
    support more mosquitoes.
    Threat Score = Surface Area * Age
    """

    ponds = {
        "A": {"size": 10, "age": 1, "description": "10 feet square, one year old."},
        "C": {"size": 30, "age": 1, "description": "30 feet square, one year old."},
        "D": {"size": 10, "age": 5, "description": "10 feet square, five years old."},
        "E": {"size": 30, "age": 5, "description": "30 feet square, five years old."}
    }

    print("Calculating the threat score for each pond...")
    print("Formula: Threat Score = (Side * Side) * Age\n")

    highest_threat = {"id": None, "score": -1}

    # Calculate and print the score for each pond
    for pid, properties in ponds.items():
        size = properties["size"]
        age = properties["age"]
        surface_area = size * size
        threat_score = surface_area * age

        print(f"Pond {pid} ({properties['description']}):")
        # Output each number in the final equation as requested
        print(f"Threat Score = ({size} * {size}) * {age} = {threat_score}\n")

        if threat_score > highest_threat["score"]:
            highest_threat["id"] = pid
            highest_threat["score"] = threat_score

    print("Conclusion:")
    print(f"Pond {highest_threat['id']} has the highest threat score of {highest_threat['score']}.")
    print("This is because it is both the largest and the oldest, providing the most ideal habitat for mosquitoes.")

solve_mosquito_threat()