def solve_hat_puzzle():
    """
    This function explains and prints the solution to the nine-person hat puzzle.
    """
    # The final distribution of hats (W for White, B for Black).
    # The list is 1-indexed for clarity, so we put a placeholder at index 0.
    hats = [None, 'W', 'B', 'B', 'W', 'W', 'B', 'B', 'W', 'B']
    total_black = 5
    total_white = 4

    print("--- The Hat Puzzle Solution ---")
    print("\nStep 1: The Initial Setup")
    print(f"There are 9 people with {total_black} Black hats and {total_white} White hats.")
    print("The final distribution of hats is:")
    config_str = " ".join(hats[1:])
    print(f"P1 P2 P3 P4 P5 P6 P7 P8 P9")
    print(f"{config_str}\n")

    print("Step 2: Information from Round 1 and 2 ('No')")
    print("After two rounds of 'No', a crucial piece of information becomes public knowledge.")
    print("The logic implies that everyone knows Person 2, 3, 6, and 7 must have Black hats.")
    
    known_black_hats = 4
    remaining_black_hats = total_black - known_black_hats
    
    print(f"\nThere are {total_black} total Black hats.")
    print(f"After Round 2, {known_black_hats} of them are publicly identified.")
    print(f"Final Equation Part 1: {total_black} - {known_black_hats} = {remaining_black_hats}")
    print(f"This means there is only {remaining_black_hats} remaining Black hat among Persons 1, 4, 5, 8, and 9.")
    print("The other four in this group must have White hats.\n")

    print("Step 3: The Deduction in Round 3 ('Yes')")
    
    yes_sayers = []
    
    # Person 4's Logic
    p4_sees_p9_hat = hats[9]
    if p4_sees_p9_hat == 'B':
        print("Person 4's logic:")
        print(" > 'I know the last Black hat is on P1, P4, P5, P8, or P9.'")
        print(" > 'I can see the hats of P1, P8, and P9.'")
        print(f" > 'I see that P9 has a {p4_sees_p9_hat} hat. That's the last Black hat!'")
        print(" > 'Therefore, my own hat must be White.'")
        yes_sayers.append(4)

    # Person 5's Logic
    p5_sees_p9_hat = hats[9]
    if p5_sees_p9_hat == 'B':
        print("\nPerson 5's logic:")
        print(" > 'My logic is the same as Person 4's.'")
        print(" > 'I also see the last Black hat on P9.'")
        print(" > 'Therefore, my own hat must also be White.'")
        yes_sayers.append(5)

    print("\n--- Final Answer ---")
    print(f"Number of people who replied 'Yes': {len(yes_sayers)}")
    print("The distribution of hats is:")
    
    # Print the final configuration clearly
    for i in range(1, 10):
        print(f"  Person {i}: {hats[i]}")

solve_hat_puzzle()