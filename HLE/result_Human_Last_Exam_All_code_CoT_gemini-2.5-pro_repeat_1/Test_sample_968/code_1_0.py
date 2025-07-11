import sys

def find_correct_arabesques():
    """
    This function identifies the Vaganova arabesques where the forward arm
    is on the opposite side of the lifted leg.
    """
    # Define the Vaganova arabesque positions.
    # The 'forward_arm_relation' key describes the forward arm's position
    # relative to the lifted leg.
    arabesques_info = {
        "First": {
            "forward_arm_relation": "opposite" # Arm of supporting leg is forward
        },
        "Second": {
            "forward_arm_relation": "same" # Arm of lifted leg is forward
        },
        "Third": {
            "forward_arm_relation": "opposite" # Arm of supporting leg is forward
        },
        "Fourth": {
            "forward_arm_relation": "same" # Arm of lifted leg is forward
        }
    }

    matching_arabesques = []
    print("Analyzing Vaganova arabesque positions:")
    for name, properties in arabesques_info.items():
        if properties["forward_arm_relation"] == "opposite":
            matching_arabesques.append(name)
            print(f"- {name} Arabesque: The forward arm is on the OPPOSITE side of the lifted leg.")
        else:
            print(f"- {name} Arabesque: The forward arm is on the SAME side as the lifted leg.")

    print("\nConclusion:")
    print(f"The two types of arabesque that fit the description are {matching_arabesques[0]} and {matching_arabesques[1]}.")
    print("This corresponds to answer choice C.")


# Execute the function to find the answer.
if __name__ == "__main__":
    find_correct_arabesques()