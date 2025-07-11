# Step 1: Define the number of curves that remained unchanged, based on the poem.
stayed_curves = 90

# Step 2: Determine the fraction of total curves that remained unchanged.
# The total is 1 (or 8/8). The disturbed fractions are 3/8 and 1/4 (or 2/8).
# Fraction that stayed = 1 - (3/8) - (1/4)
fraction_stayed = 1 - (3/8) - (2/8) # which is 3/8

# Step 3: Calculate the total number of original curves (T).
# If 3/8 of the total is 90, we can find the full total.
total_original_curves = stayed_curves / fraction_stayed

# Step 4: Calculate the number of curves in each disturbed group.
# Curves that "lost shape serene" were 3/8 of the total.
lost_shape_curves = total_original_curves * (3/8)

# Curves that "found new paths" were 1/4 of the total.
new_path_curves = total_original_curves * (1/4)

# Step 5: The total number of curves the master must draw is the sum of all disturbed curves.
total_curves_to_draw = lost_shape_curves + new_path_curves

# Step 6: Display the final equation and the answer, ensuring all numbers are integers.
print("The master must restore the curves that 'lost shape' and those that 'found new paths'.")
print("The final calculation is:")
print(f"{int(lost_shape_curves)} + {int(new_path_curves)} = {int(total_curves_to_draw)}")

print(f"\nTotal curves the master must draw: {int(total_curves_to_draw)}")
<<<150>>>