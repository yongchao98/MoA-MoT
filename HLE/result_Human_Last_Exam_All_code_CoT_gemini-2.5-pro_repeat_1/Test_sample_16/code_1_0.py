def find_minimum_prototypes(C):
    """
    Calculates the minimum number of prototypes required to guarantee that a 
    distance-weighted soft-label kNN classifier will correctly classify C class centroids.

    Args:
        C (int): The number of disjoint, contiguous, unimodal classes.

    Returns:
        int: The minimum number of prototypes (M) required.
    """
    # Based on logical deduction, the minimum number of prototypes (M) must be equal
    # to the number of classes (C).
    #
    # 1. M >= C: If M < C, the classifier can produce at most M unique class labels.
    #    This makes it impossible to guarantee a unique and correct classification for all C centroids.
    #
    # 2. M <= C: A configuration with M = C prototypes is sufficient. This can be achieved
    #    by placing one prototype exactly on each of the C class centroids.
    #
    # Therefore, the minimum number of prototypes required is exactly C.
    M = C
    return M

# Let's consider an example where there are 12 classes.
num_classes = 12

# Calculate the minimum number of prototypes required.
min_prototypes = find_minimum_prototypes(num_classes)

# The problem asks to output the numbers in the final equation.
# The final equation is: Minimum Prototypes = Number of Classes
print(f"For C = {num_classes} classes, the minimum number of prototypes required is M = {min_prototypes}.")
print("The final equation is expressed as:")
print(f"{min_prototypes} = {num_classes}")