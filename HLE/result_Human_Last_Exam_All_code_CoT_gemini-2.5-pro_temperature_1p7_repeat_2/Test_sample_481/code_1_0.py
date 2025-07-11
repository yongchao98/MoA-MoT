# Define the numerical imprints for each capsule
imprint_10mg = "3105"
imprint_20mg = "3106"

# Use set operations to find the unique digits in each imprint string.
# A set contains only unique elements.
set_10mg = set(imprint_10mg)
set_20mg = set(imprint_20mg)

# The difference() method returns a set containing the difference between two sets.
# We find the digit that is in the 10mg imprint but not in the 20mg one.
diff_1 = set_10mg.difference(set_20mg)

# We do the same for the 20mg imprint.
diff_2 = set_20mg.difference(set_10mg)

# The result of the difference operation is a set. Since we expect only one differing
# number in each, we use pop() to get the single element out of the set.
diff_num_10mg = diff_1.pop()
diff_num_20mg = diff_2.pop()

print("The imprint on the 10mg capsule contains the number: " + diff_num_10mg)
print("The imprint on the 20mg capsule contains the number: " + diff_num_20mg)
print("\nTherefore, the different numbers are:")
# As requested, outputting each number that makes up the final difference 'equation' (5 vs 6)
print(diff_num_10mg)
print(diff_num_20mg)
