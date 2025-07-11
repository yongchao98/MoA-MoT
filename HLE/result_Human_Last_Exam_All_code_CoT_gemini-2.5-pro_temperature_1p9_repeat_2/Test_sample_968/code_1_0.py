def solve_vaganova_arabesque():
    """
    Analyzes the Vaganova arabesque positions to find which ones have the
    forward arm on the opposite side of the lifted leg.
    """
    # Define the rule for each arabesque position.
    # The value describes the forward arm's side relative to the lifted leg.
    arabesque_rules = {
        "First": "opposite",
        "Second": "same",
        "Third": "opposite",
        "Fourth": "same"
    }

    print("Analyzing the four Vaganova arabesque positions:")
    print("-" * 60)

    matches = []
    for position, rule in arabesque_rules.items():
        if rule == "opposite":
            result = "MATCHES the criterion."
            matches.append(position)
        else:
            result = "Does NOT match."

        print(f"- {position} Arabesque: The forward arm is on the '{rule}' side as the lifted leg. -> {result}")

    print("-" * 60)
    # The final equation is a combination of the matching arabesques.
    if len(matches) == 2:
        print(f"The two types that satisfy the condition are the First and Third.")
        print(f"Final logical combination: {matches[0]} + {matches[1]}")
    else:
        print("Could not determine the two matching arabesques.")

solve_vaganova_arabesque()