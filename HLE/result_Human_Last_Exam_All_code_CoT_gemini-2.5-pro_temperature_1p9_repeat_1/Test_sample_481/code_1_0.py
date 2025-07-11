# The numerical imprint on Teva's 10mg fluoxetine capsule is 3105.
# The numerical imprint on Teva's 20mg fluoxetine capsule is 3106.
# Note: The text "TEVA" is also imprinted on the capsules but is excluded per the prompt.

imprint_10mg = "3105"
imprint_20mg = "3106"

# Find the digits that are not shared between the two imprints
# by converting the strings to sets and finding the symmetric difference.
diff_set = set(imprint_10mg).symmetric_difference(set(imprint_20mg))

# Convert the set of differing digits to a sorted list for consistent output.
different_digits = sorted(list(diff_set))

# Print the full numbers and then the digits that differ.
print(f"The number imprinted on the 10mg capsule is: {imprint_10mg}")
print(f"The number imprinted on the 20mg capsule is: {imprint_20mg}")
print(f"\nThe digits that are different between these two numbers are: {', '.join(different_digits)}")