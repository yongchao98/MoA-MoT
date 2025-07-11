import itertools

def solve_letter_arrangement():
    """
    Calculates the number of valid arrangements for the letters L, N, S, W
    based on a specific connection rule.
    """
    # Step 1: Define the tip types (Left, Right) for each letter.
    # 'T' for Top, 'B' for Bottom.
    # For example, L starts at the Top and ends at the Bottom.
    tip_types = {
        'L': ('T', 'B'),
        'N': ('B', 'T'),
        'S': ('T', 'B'),
        'W': ('T', 'T')
    }

    # Step 2: Define the connection rule based on the problem's examples.
    # "S can connect to W" -> Right(S)='B', Left(W)='T'. So B can connect to T.
    # "W cannot connect to S" -> Right(W)='T', Left(S)='T'. So T cannot connect to T.
    # The rule is that the tip types must be different.
    def can_connect(letter1, letter2):
        right_tip_l1 = tip_types[letter1][1]
        left_tip_l2 = tip_types[letter2][0]
        return right_tip_l1 != left_tip_l2

    letters = ['L', 'N', 'S', 'W']
    all_arrangements = list(itertools.permutations(letters))
    
    valid_arrangements = []
    
    # Step 3: Check all possible arrangements (permutations).
    for arrangement in all_arrangements:
        is_valid = True
        # Check connections between adjacent letters in the arrangement.
        if not can_connect(arrangement[0], arrangement[1]):
            is_valid = False
        if not can_connect(arrangement[1], arrangement[2]):
            is_valid = False
        if not can_connect(arrangement[2], arrangement[3]):
            is_valid = False
        
        if is_valid:
            valid_arrangements.append("".join(arrangement))

    # Step 4: Output the results.
    if valid_arrangements:
        print("Found the following valid arrangements:")
        for arr in valid_arrangements:
            # S -> L -> W -> N
            # For each connection, let's show the tip types.
            # S(B) -> L(T), L(B) -> W(T), W(T) -> N(B)
            # This follows the format of Number1 -> Number2 -> Number3 = Result.
            # Here we print the letters of the valid arrangement.
            print(f"{arr[0]} -> {arr[1]} -> {arr[2]} -> {arr[3]}")
    else:
        print("No valid arrangements found.")

    print("\nFinal equation:")
    # "Remember in the final code you still need to output each number in the final equation!"
    # This instruction can be interpreted as showing the sum of valid arrangements found.
    # Since we found 2 arrangements, the equation is 1 + 1 = 2.
    if len(valid_arrangements) > 0:
        equation_str = " + ".join(["1"] * len(valid_arrangements))
        print(f"{equation_str} = {len(valid_arrangements)}")
    else:
        print("0")
    
    print("\nTotal number of ways:")
    print(len(valid_arrangements))

solve_letter_arrangement()