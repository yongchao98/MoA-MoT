import itertools

# Step 1: Define the connection properties for each letter.
# We model the tips as being 'Top' or 'Bottom'.
# letter_tips[letter] = (Left_Tip, Right_Tip)
letter_tips = {
    'L': ('Top', 'Bottom'),
    'N': ('Bottom', 'Top'),
    'S': ('Top', 'Bottom'),
    'W': ('Top', 'Top')
}

# Step 2: Define the connection rule.
# Based on the examples S->W (valid) and W->S (invalid), we deduce the rule is that
# the tip levels must be different for a connection to be valid.
# S(right=Bottom) -> W(left=Top): valid, as Bottom != Top.
# W(right=Top)   -> S(left=Top): invalid, as Top == Top.
def is_valid_connection(letter1, letter2):
    """Checks if the right tip of letter1 can connect to the left tip of letter2."""
    right_tip_1 = letter_tips[letter1][1]
    left_tip_2 = letter_tips[letter2][0]
    return right_tip_1 != left_tip_2

# Step 3: Iterate through all unique arrangements of the four letters.
letters = ['L', 'N', 'S', 'W']
all_permutations = list(itertools.permutations(letters))

valid_arrangement_count = 0
calculation_numbers = []

for p in all_permutations:
    # Unpack the arrangement
    p1, p2, p3, p4 = p
    
    # Check the three connections in the sequence
    if (is_valid_connection(p1, p2) and
        is_valid_connection(p2, p3) and
        is_valid_connection(p3, p4)):
        
        # This is a valid arrangement
        valid_arrangement_count += 1
        # Add a '1' to our list for each valid arrangement found
        calculation_numbers.append(1)

# Step 4: Output the result as an equation, as requested.
# If two arrangements are found, this will print "1 + 1 = 2".
if valid_arrangement_count > 0:
    equation_str = " + ".join(str(num) for num in calculation_numbers)
    print(f"{equation_str} = {valid_arrangement_count}")
else:
    # If no arrangements are found, the answer is 0.
    print(0)