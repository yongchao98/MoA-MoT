# The analysis identified three triplets where all components are directly related:
# 1. Aristotelian, Charadrius vociferus - wing, Recurvirostra avosetta - wing
# 5. Müllerian, Heliconiini - color, Melipotini - color
# 8. Aposematism, Danaus plexipus - wing, Cycnia tenera - tymbal

# Storing the numbers of the correct triplets
correct_triplet_numbers = [1, 5, 8]

# Sorting the numbers in ascending order
correct_triplet_numbers.sort()

# Printing the final numbers as a comma-separated string
output = ", ".join(map(str, correct_triplet_numbers))
print("The numbers of the correct triplets in ascending order are:")
print(output)