# Initial parameters
initial_rope_length = 60
shrink_factor = 0.75

# Step 1: Define the first cut.
# We choose a cut that is a multiple of 4 and more than half the total length.
# Let's test the scenario where the first cut creates a 40cm piece.
longer_portion_first_cut = 40
shorter_portion_first_cut = initial_rope_length - longer_portion_first_cut

# Step 2: Calculate the length of the longer portion after it shrinks.
# This is the answer to the problem's question.
longer_portion_after_shrink = longer_portion_first_cut * shrink_factor

# Step 3: To prove this is a valid solution, we show how the second cut results in a 15cm piece.
# The two pieces we have are longer_portion_after_shrink (30cm) and shorter_portion_first_cut (20cm).
# We can cut the 30cm piece to get a 15cm piece.
# If the shrunken piece of the second cut is 15cm, its original length must have been 15 / 0.75 = 20cm.
longer_portion_second_cut = 20

# The piece we cut (30cm) is split into the longer_portion_second_cut (20cm) and a remainder.
shorter_portion_second_cut = longer_portion_after_shrink - longer_portion_second_cut

# The longer portion of the second cut (20cm) shrinks.
final_piece_from_shrink = longer_portion_second_cut * shrink_factor

# Print the step-by-step reasoning and the final calculation.
print("Problem: Find the length of the longer portion after the first cut.")
print(f"1. A 60cm rope is cut into a longer piece of {longer_portion_first_cut}cm and a shorter piece of {shorter_portion_first_cut}cm.")
print(f"2. The longer piece ({longer_portion_first_cut}cm) shrinks by 25%.")
print(f"3. The calculation for its new length is: {longer_portion_first_cut} * {shrink_factor} = {int(longer_portion_after_shrink)}")
print(f"The length of the longer portion after the first cut is {int(longer_portion_after_shrink)}cm.")
print("\nVerification:")
print(f"To verify, we take the {int(longer_portion_after_shrink)}cm piece and cut it into {longer_portion_second_cut}cm and {int(shorter_portion_second_cut)}cm.")
print(f"The new longer piece ({longer_portion_second_cut}cm) shrinks by 25%, resulting in a {int(final_piece_from_shrink)}cm piece.")
print("This confirms it's possible to get a 15cm rope with this initial cut.")
