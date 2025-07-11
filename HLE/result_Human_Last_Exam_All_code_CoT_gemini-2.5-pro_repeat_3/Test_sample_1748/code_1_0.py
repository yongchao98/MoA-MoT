def solve_rope_puzzle():
    """
    This function solves the magical rope puzzle by finding a valid sequence of cuts.
    """
    initial_rope_length = 60
    target_piece = 15
    shrink_rate = 0.75
    
    # We are looking for p1_long, the length of the longer piece from the first cut.
    # It must be greater than half the rope's length (30) and be a multiple of 4
    # for its shrunken length to be an integer.
    # We will search for the smallest possible solution.
    for p1_long in range(32, initial_rope_length, 4):
        
        # From our derived equation: p2_long = 45 - (p1_long / 4)
        # where p2_long is the longer piece from the second cut.
        p2_long = 45 - (p1_long / 4)
        
        # Check if this p2_long is a valid solution.
        # It must be an integer, greater than the target piece (15),
        # and a multiple of 4 for its own shrunken length to be an integer.
        is_integer = (p2_long == int(p2_long))
        is_longer = (p2_long > target_piece)
        is_multiple_of_4 = (p2_long % 4 == 0)
        
        if is_integer and is_longer and is_multiple_of_4:
            p1_short = initial_rope_length - p1_long
            p1_long_shrunk = int(p1_long * shrink_rate)
            
            # This is the length of the longer portion *after* the first cut's shrink.
            final_answer = p1_long_shrunk

            print("Found a valid scenario:")
            print(f"1. First cut: The {initial_rope_length}cm rope is cut into a {p1_short}cm piece and a {p1_long}cm piece.")
            print(f"2. Shrinking: The longer piece ({p1_long}cm) shrinks by 25% to become {p1_long_shrunk}cm.")
            
            new_total_length = p1_short + p1_long_shrunk
            print(f"3. New rope length for the second cut: {p1_short} + {p1_long_shrunk} = {new_total_length}cm.")

            p2_long = int(p2_long)
            p2_short = new_total_length - p2_long
            print(f"4. Second cut: The {new_total_length}cm rope is cut into a {p2_short}cm piece and a {p2_long}cm piece.")

            print("\nThe final equation derived from the problem constraints is:")
            print(f"45 - (p1_long / 4) = p2_long")
            print(f"The numbers that solve this equation are p1_long = {p1_long} and p2_long = {p2_long}.")
            print(f"Verification: 45 - ({p1_long} / 4) = 45 - {int(p1_long/4)} = {45 - int(p1_long/4)}, which equals {p2_long}.")

            print(f"\nThe question asks for the length of the longer portion after the first cut.")
            print(f"This is the length of the {p1_long}cm piece after it shrinks: {p1_long} * {shrink_rate} = {final_answer}")
            
            # Stop after finding the first solution
            return final_answer

# Run the solver and capture the result
result = solve_rope_puzzle()
# The final answer is printed in the required format below
print(f"\n<<<__{result}__>>>")
