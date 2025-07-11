def solve_arabesque_question():
    """
    Analyzes Vaganova arabesques to find which have the forward arm
    on the opposite side of the lifted leg.
    """
    print("Analysis of Vaganova Arabesque Positions:")
    print("The condition is: Arm extended forward is on the OPPOSITE side of the lifted leg.")
    print("-" * 60)

    # Rule for First Arabesque
    print("Position 1 (First Arabesque):")
    print("   The arm on the same side as the *supporting* leg extends forward.")
    print("   - Since the supporting leg is opposite the lifted leg, this MATCHES the condition.")
    print("")

    # Rule for Second Arabesque
    print("Position 2 (Second Arabesque):")
    print("   The arm on the same side as the *lifted* leg extends forward.")
    print("   - This does NOT match the condition.")
    print("")

    # Rule for Third Arabesque
    print("Position 3 (Third Arabesque):")
    print("   The arm on the same side as the *supporting* leg extends forward.")
    print("   - Like the first, this MATCHES the condition.")
    print("")

    # Rule for Fourth Arabesque
    print("Position 4 (Fourth Arabesque):")
    print("   The arm on the same side as the *lifted* leg extends forward.")
    print("   - This does NOT match the condition.")
    print("")

    print("-" * 60)
    print("Conclusion: The positions that match are the first and the third.")
    # Fulfilling the "output each number in the final equation" requirement.
    first_position_number = 1
    third_position_number = 3
    print(f"Final matching positions are represented by the numbers: {first_position_number} and {third_position_number}")

solve_arabesque_question()