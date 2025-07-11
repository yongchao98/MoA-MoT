import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters L, N, S, W can be arranged
    based on a specific connection rule for their tips.
    """

    # Step 1: Define the properties of each letter's tips.
    # The tuple represents (leftmost_tip_type, rightmost_tip_type).
    letter_tips = {
        'L': ('TOP', 'BOTTOM'),
        'N': ('BOTTOM', 'TOP'),
        'S': ('TOP', 'BOTTOM'),
        'W': ('TOP', 'TOP')
    }
    
    letters = list(letter_tips.keys())

    # Step 2: Define the connection rule function based on our deduction.
    # A connection from letter1 to letter2 is valid if the rightmost tip of letter1
    # is of a different type than the leftmost tip of letter2.
    def can_connect(letter1, letter2):
        right_tip1 = letter_tips[letter1][1]
        left_tip2 = letter_tips[letter2][0]
        return right_tip1 != left_tip2

    # Step 3: Iterate through all permutations and count the valid ones.
    valid_arrangements = []
    all_permutations = itertools.permutations(letters)

    for p in all_permutations:
        is_valid = True
        # Check connections for adjacent letters in the permutation.
        for i in range(len(p) - 1):
            if not can_connect(p[i], p[i+1]):
                is_valid = False
                break
        
        if is_valid:
            valid_arrangements.append("".join(p))
    
    # Step 4: Output the result as an equation.
    count = len(valid_arrangements)
    
    if count > 0:
        # Create a string like "1 + 1 = 2"
        sum_parts = ["1"] * count
        equation_str = f"{' + '.join(sum_parts)} = {count}"
        print(equation_str)
    else:
        # If no arrangements are found.
        print("0")

solve_letter_arrangement()