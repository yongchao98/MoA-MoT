def find_matching_arabesques():
    """
    Identifies Vaganova arabesques based on arm and leg positions.
    """
    # In Vaganova terminology:
    # 'opposite' means the forward arm is on the opposite side of the lifted leg.
    # 'same' means the forward arm is on the same side as the lifted leg.
    vaganova_arabesques = {
        1: "opposite",
        2: "same",
        3: "opposite",
        4: "same"
    }

    matching_numbers = []
    for number, arm_position in vaganova_arabesques.items():
        if arm_position == "opposite":
            matching_numbers.append(number)

    # Ensure the numbers are sorted for a clear output
    matching_numbers.sort()

    print("Analyzing the four Vaganova arabesques:")
    print("- In the First Arabesque, the arm extended forward is on the opposite side of the lifted leg.")
    print("- In the Second Arabesque, the arm extended forward is on the same side as the lifted leg.")
    print("- In the Third Arabesque, the arm extended forward is on the opposite side of the lifted leg.")
    print("- In the Fourth Arabesque, the arm extended forward is on the same side as the lifted leg.")
    print("\nBased on this, the two types of arabesque matching the criteria are:")
    
    # The prompt requires outputting each number in the final equation/statement.
    print(f"The {matching_numbers[0]}st and {matching_numbers[1]}rd.")

if __name__ == "__main__":
    find_matching_arabesques()
<<<C>>>