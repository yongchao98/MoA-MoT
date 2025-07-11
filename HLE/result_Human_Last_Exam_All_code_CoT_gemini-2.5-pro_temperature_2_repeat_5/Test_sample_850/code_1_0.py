# The user wants to find crystal classes that satisfy three conditions:
# 1. Achiral
# 2. Non-polar
# 3. Have the correct symmetry for optical activity

# The "correct symmetry for optical activity" is chirality.
# By definition, a chiral object is not superimposable on its mirror image.
# This means chiral crystal classes lack improper rotation axes (like inversion centers or mirror planes).

# An achiral object IS superimposable on its mirror image, meaning it MUST have an
# improper rotation axis (e.g., a mirror plane or inversion center).

# Therefore, the conditions 'achiral' and 'optically active (chiral)' are mutually exclusive.
# This script will demonstrate that there is no solution by finding the intersection of the relevant sets of crystal classes.

# Define the sets of crystal classes based on their properties. We use Hermann-Mauguin notation.
optically_active_classes = {
    '1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'
}

# The set of all 32 crystal classes
all_32_classes = {
    '1', 'overline{1}', '2', 'm', '2/m', '222', 'mm2', 'mmm', '3', 'overline{3}',
    '32', '3m', 'overline{3}m', '4', 'overline{4}', '4/m', '422', '4mm', 'overline{4}2m',
    '4/mmm', '6', 'overline{6}', '6/m', '622', '6mm', 'overline{6}m2', '6/mmm',
    '23', 'm-3', '432', 'overline{4}3m', 'm-3m'
}

# Achiral classes are all classes that are NOT optically active (chiral).
achiral_classes = all_32_classes - optically_active_classes

# The 10 polar classes exhibit a unique polar direction.
polar_classes = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}

# Non-polar classes are all classes that are NOT polar.
non_polar_classes = all_32_classes - polar_classes

# Now we find the intersection of the three requested sets.
# The "equation" is: (Achiral) & (Non-Polar) & (Optically Active)
result_set = achiral_classes.intersection(non_polar_classes, optically_active_classes)

# Print the explanation and the result.
print("--- The Problem ---")
print("To find crystal classes that are: Achiral AND Non-Polar AND Optically Active.")
print("\n--- The Physics ---")
print("A crystal must be CHIRAL to be Optically Active. By definition, a CHIRAL crystal cannot be ACHIRAL.")
print("The request contains a contradiction. Let's prove this by showing the sets and their intersection.")

# Print the members ("numbers") of each set in the "equation"
print("\n--- Equation Sets ---")
print("\nSet 1: Achiral Classes")
print(' '.join(sorted(list(achiral_classes))))

print("\nSet 2: Non-Polar Classes")
print(' '.join(sorted(list(non_polar_classes))))

print("\nSet 3: Optically Active (Chiral) Classes")
print(' '.join(sorted(list(optically_active_classes))))


# Print the final result of the intersection.
print("\n--- Final Result ---")
if not result_set:
    print("The intersection of these three sets is EMPTY.")
    print("There are NO crystal classes that are achiral, non-polar, and optically active.")
else:
    # This code path is unreachable due to the logic of the physics, but is included for completeness.
    print("The following classes meet the criteria: " + ', '.join(sorted(list(result_set))))
