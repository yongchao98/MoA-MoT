def find_minimum_prototypes(C, N=None, D=None):
    """
    Determines the minimum number of prototypes required to guarantee that a
    distance-weighted soft-label kNN classifier can correctly classify the
    centroids of C distinct classes.

    The logic proceeds in two steps:
    1. Proving that C-1 prototypes are insufficient.
    2. Proving that C prototypes are sufficient.

    Args:
        C (int): The number of disjoint classes.
        N (int, optional): The total number of datapoints. This parameter is not
                         needed for the theoretical solution but is included for
                         completeness based on the problem description.
        D (int, optional): The number of dimensions of the manifold. This is also
                         not needed for the theoretical solution.
    """

    # --- Step 1: Prove Necessity (M >= C) ---
    # We must show that any number of prototypes M < C is not enough to guarantee
    # success for *all* possible arrangements of the class centroids.
    #
    # Consider the case with M = C-1 prototypes. We can create a worst-case
    # scenario for these M prototypes.
    # Let's place the first C-1 centroids (µ_1, ..., µ_{C-1}) exactly where the M
    # prototypes (p_1, ..., p_M) are. Let p_i = µ_i.
    # The soft label for prototype p_i is based on its proximity to all class
    # centroids. Since d(p_i, µ_i) = 0, its proximity is infinite, and its soft
    # label effectively becomes a "hard label" for class i: s_i = (0,..,1,..0).
    #
    # Now, place the last centroid, µ_C, extremely far away from all C-1 prototypes.
    #
    # To classify µ_C, the kNN classifier finds its k nearest prototypes. All of
    # these will be from the set {p_1, ..., p_{C-1}}. The soft labels for these
    # prototypes are s_1, ..., s_{C-1}, which represent classes 1 to C-1.
    # None of these soft labels have a non-zero value for the C-th class.
    # The final predicted soft label for µ_C will be a weighted average of vectors
    # that are all zero in the C-th position. Thus, the prediction for µ_C
    # will also have a zero in the C-th position, making correct classification
    # as class C impossible.
    #
    # This proves that M = C-1 prototypes is not enough. Therefore, we need at least C prototypes.

    # --- Step 2: Prove Sufficiency (M = C) ---
    # We must show that M = C prototypes are enough to guarantee success.
    # We can do this with a constructive placement strategy.
    #
    # Place each of the C prototypes exactly at the location of one of the C
    # class centroids. That is, for i=1 to C, set prototype p_i = µ_i.
    #
    # As established before, the soft label for prototype p_i will be a hard
    # label for class i (s_i has a 1 at position i and 0s elsewhere).
    #
    # Now, let's classify any centroid µ_j. The classifier measures the distance from
    # µ_j to all C prototypes. The prototype p_j is at distance 0, while all other
    # prototypes p_i (where i != j) are at a distance d(µ_j, µ_i) > 0.
    #
    # The kNN algorithm will identify p_j as the nearest prototype. In the case of
    # an exact match, the classifier assigns the label of that point. The label of
    # prototype p_j corresponds to class j.
    # This means µ_j will be correctly classified. Since this works for any j,
    # all C centroids will be classified correctly.

    # --- Conclusion ---
    # Since C-1 prototypes are not sufficient, and C prototypes are sufficient,
    # the minimum number required is C. The number of datapoints N and the
    # dimensionality D do not affect this worst-case theoretical guarantee.

    minimum_number_of_prototypes = C
    
    # Final equation: Minimum Prototypes = C
    print(f"The number of classes is C = {C}.")
    print("The final equation determining the minimum number of prototypes is:")
    print(f"Minimum Prototypes = {minimum_number_of_prototypes}")


# --- Execution ---
# You can define the number of classes C and run the function.
# For example, if you have 5 classes:
number_of_classes = 5
find_minimum_prototypes(C=number_of_classes)