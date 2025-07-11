# Set initial parameters based on the problem
initial_rope_length = 60
shrink_factor = 0.75
target_piece = 15

# --- First Cut ---
# Based on our analysis, the most elegant solution involves cutting the rope into 1/3 and 2/3.
first_cut_short = initial_rope_length / 3
first_cut_long = initial_rope_length * 2 / 3

l1 = int(first_cut_long)
s1 = int(first_cut_short)

print(f"1. Start with a rope of {initial_rope_length} cm.")
print(f"2. Make the first cut, creating pieces of {s1} cm and {l1} cm.")

# --- First Shrinkage ---
l1_shrunk = int(l1 * shrink_factor)
print(f"3. The longer piece ({l1} cm) shrinks by 25% to {l1_shrunk} cm.")
rope_for_second_cut = l1_shrunk

# --- Second Cut ---
# To get our target piece, we need a pre-shrink length of 20cm (since 20 * 0.75 = 15)
# We can get a 20cm piece by cutting the 30cm piece.
second_cut_long = int(target_piece / shrink_factor)
second_cut_short = rope_for_second_cut - second_cut_long

l2 = second_cut_long
s2 = second_cut_short

print(f"4. Take the {rope_for_second_cut} cm piece and cut it into {s2} cm and {l2} cm.")

# --- Second Shrinkage ---
l2_shrunk = int(l2 * shrink_factor)
print(f"5. The longer piece from this cut ({l2} cm) shrinks by 25% to {l2_shrunk} cm.")

print(f"\nWe have successfully obtained a {l2_shrunk} cm piece in two cuts.")
print(f"The length of the longer portion after the first cut was {l1} cm.")

# Final answer: l1
# Final Equation: ((60 * 2/3) * 3/4) is the piece used for the second cut
# ((30) * 2/3) * 3/4 = 15
print(f"\nFinal Equation: (({rope_for_second_cut}) * {l2}/{rope_for_second_cut}) * {shrink_factor} = {target_piece}")
