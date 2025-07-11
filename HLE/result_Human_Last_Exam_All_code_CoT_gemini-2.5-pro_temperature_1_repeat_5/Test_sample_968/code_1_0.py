def find_matching_arabesques():
    """
    This script analyzes the Vaganova arabesque positions to find which ones
    have the forward arm on the opposite side of the lifted leg.
    """

    # We define the rules for each arabesque.
    # The key is the arabesque number, and the value describes the forward arm's position
    # relative to the legs.
    # "Opposite the lifted leg" is the same as "same side as the supporting leg".
    arabesque_rules = {
        1: "Forward arm is on the same side as the supporting leg (opposite the lifted leg).",
        2: "Forward arm is on the same side as the lifted leg.",
        3: "The arm on the side of the supporting leg is extended forward and high (opposite the lifted leg).",
        4: "Forward arm is on the same side as the lifted leg."
    }

    matching_arabesques = []
    print("Analyzing Vaganova arabesques...")
    print("Condition: Forward arm is on the opposite side of the lifted leg.\n")

    # Iterate through the rules and check for the condition.
    for number, rule in sorted(arabesque_rules.items()):
        # We check if the rule explicitly states the arm is opposite the lifted leg
        # or on the same side as the supporting leg.
        if "opposite the lifted leg" in rule or "side of the supporting leg" in rule:
            matching_arabesques.append(number)
            print(f" - Arabesque {number}: Matches. {rule}")
        else:
            print(f" - Arabesque {number}: Does not match. {rule}")

    # Sort the results for a clean output.
    matching_arabesques.sort()
    
    number1 = matching_arabesques[0]
    number2 = matching_arabesques[1]

    print("\nConclusion:")
    print(f"The two arabesques that satisfy the condition are the First and Third.")
    print(f"The numbers in the final answer are {number1} and {number2}.")

find_matching_arabesques()
<<<C>>>