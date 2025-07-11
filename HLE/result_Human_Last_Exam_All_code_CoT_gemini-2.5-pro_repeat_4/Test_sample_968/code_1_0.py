def solve_arabesque_question():
    """
    This function determines which Vaganova arabesques have the forward arm
    on the opposite side of the lifted leg.
    """
    # In the Vaganova method:
    # - Supporting leg: The leg the dancer is standing on.
    # - Lifted leg: The leg extended behind the dancer.
    # The "opposite side as the lifted leg" is the same as the "same side as the supporting leg".

    vaganova_rules = {
        "First": "opposite",  # Forward arm is on the same side as the supporting leg.
        "Second": "same",     # Forward arm is on the same side as the lifted leg.
        "Third": "opposite",  # Both arms are forward, with the one corresponding to the supporting leg being higher.
        "Fourth": "same"      # Forward arm is on the same side as the lifted leg, with the body in a strong twist.
    }

    print("Analyzing the four Vaganova arabesques based on the position of the forward arm relative to the lifted leg:")
    
    matching_arabesques = []
    for arabesque, arm_relation in vaganova_rules.items():
        if arm_relation == "opposite":
            matching_arabesques.append(arabesque)
            print(f"- {arabesque} Arabesque: The forward arm is on the OPPOSITE side of the lifted leg. (Matches criteria)")
        else:
            print(f"- {arabesque} Arabesque: The forward arm is on the SAME side as the lifted leg.")

    print(f"\nThe two types of arabesque that fit the criteria are the {matching_arabesques[0]} and {matching_arabesques[1]}.")
    print("This corresponds to option C.")

solve_arabesque_question()