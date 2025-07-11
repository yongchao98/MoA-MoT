# The numerical imprint on a Teva 10mg fluoxetine capsule is 3107.
imprint_10mg = "3107"

# The numerical imprint on a Teva 20mg fluoxetine capsule is 3108.
imprint_20mg = "3108"

# We will find the numbers that are different between the two imprints.
diff_num_1 = ''
diff_num_2 = ''

# Iterate through the strings to find the first point of difference.
for i in range(len(imprint_10mg)):
    if imprint_10mg[i] != imprint_20mg[i]:
        diff_num_1 = imprint_10mg[i]
        diff_num_2 = imprint_20mg[i]
        break

# Display the full numbers and the differing digits.
# The numbers on the capsules are 3107 and 3108.
# The different numbers are 7 and 8.
print(f"The number imprinted on the 10mg capsule is {imprint_10mg[0]}{imprint_10mg[1]}{imprint_10mg[2]}{imprint_10mg[3]}.")
print(f"The number imprinted on the 20mg capsule is {imprint_20mg[0]}{imprint_20mg[1]}{imprint_20mg[2]}{imprint_20mg[3]}.")
print(f"The different numbers are {diff_num_1} and {diff_num_2}.")