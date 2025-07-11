def solve_mancala_puzzle():
    """
    Solves the Mancala puzzle by calculating the total number of stones
    and using the parity principle to determine impossible score differences.
    """
    # Define the initial game state
    p1_pits = [0, 2, 0, 0, 2, 0]
    p1_store = 22
    p2_pits = [1, 0, 0, 0, 0, 0]
    p2_store = 21

    # Calculate the sum of stones for each player's pits
    p1_pit_sum = sum(p1_pits)
    p2_pit_sum = sum(p2_pits)

    # Calculate the total number of stones in the game
    total_stones = p1_pit_sum + p1_store + p2_pit_sum + p2_store

    # --- Step 1: Display the calculation for the total number of stones ---
    print("Step 1: Calculate the total number of stones in the game.")
    # Create the string for the equation showing all parts
    p1_pit_str = " + ".join(map(str, p1_pits))
    p2_pit_str = " + ".join(map(str, p2_pits))
    print("The total is the sum of all stones on the board:")
    print(f"({p1_pit_str}) + {p1_store} + ({p2_pit_str}) + {p2_store} = {total_stones}")
    print("-" * 20)

    # --- Step 2: Explain the parity principle ---
    print("Step 2: Apply the parity principle.")
    print("Let the final scores be S1 and S2. The sum of the scores must be the total number of stones.")
    print(f"S1 + S2 = {total_stones}")
    print("\nThe score difference is D = |S1 - S2|.")
    print("A mathematical rule states that (S1 + S2) and (S1 - S2) must have the same parity (both even or both odd).")

    if total_stones % 2 == 0:
        print(f"\nSince the total number of stones ({total_stones}) is EVEN, the score difference D must also be EVEN.")
        required_parity = "even"
    else:
        print(f"\nSince the total number of stones ({total_stones}) is ODD, the score difference D must also be ODD.")
        required_parity = "odd"
    print("-" * 20)

    # --- Step 3: Evaluate the answer choices ---
    print("Step 3: Evaluate the given answer choices.")
    answer_choices = {
        "A": 0, "B": 1, "C": 2,
        "D": 3, "E": 4, "F": 5
    }

    unobtainable_choices = []
    for choice, difference in answer_choices.items():
        is_even = (difference % 2 == 0)
        if required_parity == "even" and not is_even:
            print(f"- Difference {difference} (Choice {choice}) is ODD. This is NOT possible.")
            unobtainable_choices.append(choice)
        elif required_parity == "odd" and is_even:
            print(f"- Difference {difference} (Choice {choice}) is EVEN. This is NOT possible.")
            unobtainable_choices.append(choice)
        else:
            print(f"- Difference {difference} (Choice {choice}) is {required_parity.upper()}. This is theoretically possible.")
    print("-" * 20)
    
    # --- Step 4: Final Conclusion ---
    print("Step 4: Formulate the final answer.")
    if len(unobtainable_choices) > 1:
        print(f"There are {len(unobtainable_choices)} unobtainable score differences in the list: {', '.join(unobtainable_choices)}.")
        print("Therefore, the correct option is G: 'More than one of the listed score differences is unobtainable'.")
    elif len(unobtainable_choices) == 1:
        print(f"The only unobtainable score difference is choice {unobtainable_choices[0]}.")
    else:
        print("All listed choices are theoretically possible.")

# Execute the function to solve the puzzle
solve_mancala_puzzle()