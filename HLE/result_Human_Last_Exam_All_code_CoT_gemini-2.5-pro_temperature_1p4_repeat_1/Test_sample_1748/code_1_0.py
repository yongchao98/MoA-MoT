import sys

def solve_rope_riddle():
    """
    Solves the magical rope riddle by working backwards from the desired outcome.
    """
    # Initial parameters from the problem statement
    initial_rope_length = 60
    target_piece_length = 15
    shrink_percentage = 25

    # --- Step 1: Determine the length needed for the second cut ---
    # The simplest way to get a 15cm piece is to cut a piece of rope in half.
    # This means we need a piece of 15 * 2 = 30cm for our second cut.
    # If we cut a 30cm piece into two 15cm pieces, they are of equal length,
    # and the "longer portion shrinks" rule does not apply.
    length_for_second_cut = target_piece_length * 2

    # --- Step 2: Determine the state after the first cut ---
    # The question asks for the length of the longer portion AFTER the first cut.
    # Based on our logic, this length must be 30cm.
    # Let's verify this by calculating the dimensions of the first cut.
    length_after_first_cut = length_for_second_cut

    # --- Step 3: Calculate the original length of the longer piece before it shrank ---
    # Let 'y' be the length of the longer piece just after cutting, before it shrinks.
    # The equation is: y * (1 - shrink_percentage / 100) = length_after_first_cut
    shrink_factor = 1 - (shrink_percentage / 100)
    # y = length_after_first_cut / shrink_factor
    original_longer_length = length_after_first_cut / shrink_factor

    # --- Step 4: Calculate the shorter piece from the first cut ---
    original_shorter_length = initial_rope_length - original_longer_length

    # --- Step 5: Print the full story and the final equation ---
    print("### The Logic Step-by-Step ###")
    print(f"1. To get a {target_piece_length}cm piece, we aim to cut a {length_for_second_cut}cm piece in half. This way, no piece is 'longer' and the shrinkage rule is bypassed.")
    print(f"2. This means after the first cut, the longer portion of the rope must shrink to become {length_for_second_cut}cm long.")
    print(f"3. The initial rope is {initial_rope_length}cm. It is cut into a longer piece (y) and a shorter piece (x).")
    print(f"4. The longer piece (y) shrinks by {shrink_percentage}%. To find its original length, we solve the equation:")

    # As requested, outputting the numbers in the final equation for clarity
    print("\n### The Final Equation ###")
    equation_part_1 = int(original_longer_length)
    equation_part_2 = shrink_factor
    equation_result = int(length_after_first_cut)

    print(f"The equation that gives the final length is:")
    print(f"{equation_part_1} * (1 - {shrink_percentage/100}) = {equation_result}")

    print("\nThe numbers in this equation are:")
    # Printing each number as requested
    print(equation_part_1)
    print(1)
    print(shrink_percentage / 100)
    print(equation_result)

    print("\n### Conclusion ###")
    print(f"The first cut divides the {initial_rope_length}cm rope into a {int(original_longer_length)}cm piece and a {int(original_shorter_length)}cm piece.")
    print(f"The longer portion ({int(original_longer_length)}cm) then shrinks, and its length becomes {equation_result}cm.")
    print("\nTherefore, the length of the longer portion after the first cut is:")
    print(equation_result)


solve_rope_riddle()