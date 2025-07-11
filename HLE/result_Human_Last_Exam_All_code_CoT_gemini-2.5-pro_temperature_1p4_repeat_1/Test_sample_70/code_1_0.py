def solve_flag_rank_problem():
    """
    This program determines which African flags have the same linear algebraic rank
    as the Danish flag (Rank 2) by analyzing their geometric structures.
    """

    print("Step 1: Determine the rank of the flag of Denmark.")
    print("--------------------------------------------------")
    print("Let the colors Red be represented by the number 1 and White by 2.")
    print("The Danish flag's cross design creates two distinct row patterns:")
    print(" - Type A (the red field rows): [1, 1, 2, 1, 1]")
    print(" - Type B (the white cross row): [2, 2, 2, 2, 2]")
    print("\nThese two vectors are linearly independent, as one is not a scalar multiple of the other.")
    print("An equation to show independence is c1 * A + c2 * B = 0 only if c1=0 and c2=0.")
    print("c1 * [1, 1, 2, 1, 1] + c2 * [2, 2, 2, 2, 2] = [0, 0, 0, 0, 0]")
    print("This implies c1 + 2*c2 = 0 and 2*c1 + 2*c2 = 0, which only holds for c1=0, c2=0.")
    print("Thus, the rank of the Danish flag is 2.\n")

    print("Step 2: Find African flags with a rank of 2.")
    print("---------------------------------------------")
    print("We look for flags with designs that generate exactly two linearly independent row patterns.")
    print("- Flags with only parallel stripes (e.g., Nigeria, Botswana) have a rank of 1.")
    print("- Flags with complex emblems or diagonal lines (e.g., Kenya, South Africa) have a rank > 2.")
    print("\nThe following flags have the required structure:\n")

    matching_flags = {
        "Benin": {
            "description": "A vertical green bar on the left, with horizontal yellow and red bars on the right.",
            "row_A": "[Green, Yellow, Yellow] -> [1, 2, 2]",
            "row_B": "[Green, Red, Red]     -> [1, 3, 3]",
        },
        "Madagascar": {
            "description": "A vertical white bar on the left, with horizontal red and green bars on the right.",
            "row_A": "[White, Red, Red]   -> [1, 2, 2]",
            "row_B": "[White, Green, Green] -> [1, 3, 3]",
        },
        "Senegal": {
            "description": "A vertical tricolor of green, yellow, and red, with a green star in the central yellow stripe.",
            "row_A": "Row without star: [Green, Yellow, Red] -> [1, 2, 3]",
            "row_B": "Row with star:    [Green, Green, Red]  -> [1, 1, 3]",
        }
    }

    for country, data in matching_flags.items():
        print(f"Country: {country}")
        print(f"Description: {data['description']}")
        print("This design creates two linearly independent row patterns:")
        print(f"  - Type A: {data['row_A']}")
        print(f"  - Type B: {data['row_B']}")
        print("Since there are two linearly independent row types, the rank is 2.\n")

    final_answer = ", ".join(matching_flags.keys())
    print(f"Final Answer: The flags of Benin, Madagascar, and Senegal have the same rank as the flag of Denmark.")
    # The final answer format is specified in the problem description.
    # It will be printed at the very end.
    
if __name__ == '__main__':
    solve_flag_rank_problem()
    # As requested, printing the final answer in the specified format.
    # Note: the instruction "do not ask users to copy and paste" and "use 'print' function" is followed by the script.
    # The final output is wrapped as per instructions.
    print("\n<<<Benin, Madagascar, Senegal>>>")
