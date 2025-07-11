def solve_arabesque_question():
    """
    This function identifies the Vaganova arabesques where the forward arm
    is on the opposite side of the lifted leg.
    """
    
    # In ballet terminology, "the arm opposite the lifted leg" is the same as
    # "the arm on the same side as the supporting leg". We will use this fact to identify the correct positions.
    vaganova_arabesques = {
        "First": {
            "number": 1,
            "forward_arm_is_same_as_supporting_leg": True
        },
        "Second": {
            "number": 2,
            "forward_arm_is_same_as_supporting_leg": False
        },
        "Third": {
            "number": 3,
            "forward_arm_is_same_as_supporting_leg": True
        },
        "Fourth": {
            "number": 4,
            "forward_arm_is_same_as_supporting_leg": False
        }
    }

    print("Analyzing which Vaganova arabesques have the forward arm on the opposite side of the lifted leg...")
    
    correct_names = []
    correct_numbers = []

    for name, properties in vaganova_arabesques.items():
        if properties["forward_arm_is_same_as_supporting_leg"]:
            correct_names.append(name)
            correct_numbers.append(str(properties["number"]))

    print("\nThe arabesques that fit the description are:")
    for name in correct_names:
        print(f"- {name} Arabesque")

    # This fulfills the instruction: "you still need to output each number in the final equation!"
    print("\nThe equation of the correct arabesque numbers is:")
    print(f"Correct Combination = {' + '.join(correct_numbers)}")

solve_arabesque_question()
<<<C>>>