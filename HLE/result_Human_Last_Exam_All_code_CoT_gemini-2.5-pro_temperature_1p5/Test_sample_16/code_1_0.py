import textwrap

def explain_minimum_prototypes(C):
    """
    Explains the logic for determining the minimum number of prototypes
    required for a distance-weighted kNN classifier to correctly classify
    all C class centroids.

    Args:
        C (int): The number of classes. Must be a positive integer.
    """
    if not isinstance(C, int) or C < 1:
        print("Error: The number of classes (C) must be a positive integer.")
        return

    # --- Introduction ---
    print(f"Analysis for C = {C} classes")
    print("=" * 60)
    print("Question: What is the minimum number of prototypes required to guarantee")
    print("that a distance-weighted soft-label kNN classifier will correctly")
    print("classify each of the C class centroids?")
    print("-" * 60)

    # --- Part 1: Sufficiency (C prototypes are enough) ---
    print("\nPart 1: Proving C prototypes are SUFFICIENT")
    explanation_part1 = f"""
    We can construct a set of C prototypes that guarantees success. The
    strategy is as follows:
    1. For each class `i` (from 1 to {C}), create one prototype, `P_i`.
    2. Place prototype `P_i` at the exact location of its corresponding
       class centroid, `CC_i`.
    3. Assign a 'pure' soft label to `P_i`. This label is a vector
       with a 1 in the i-th position and 0s everywhere else.
       (e.g., for C={C}, the label for `P_1` is [1, 0, ..., 0]).

    This strategy works because when we classify a centroid `CC_j`, its
    nearest prototype is `P_j` (at a distance of 0). In a
    distance-weighted kNN (where weight is 1/distance), the weight for
    `P_j` becomes infinite, while all other prototype weights are finite.
    This infinite weight means the predicted label for `CC_j` is entirely
    determined by `P_j`'s label, which correctly identifies class `j`. This
    holds true for any `k >= 1`.
    """
    print(textwrap.dedent(explanation_part1))
    print("--> Conclusion for Part 1: C prototypes are sufficient.\n")


    # --- Part 2: Necessity (Fewer than C are not enough) ---
    print("Part 2: Proving C prototypes are NECESSARY")
    explanation_part2 = f"""
    We can prove that any number of prototypes less than C is
    insufficient. Let's assume we use M prototypes, where M < C (for
    this example, M = {C - 1} prototypes). To provide a guarantee, our
    setup must work for any k; if it fails for even one k, the
    guarantee is void. Let's test the simplest case, k=1.

    With k=1, the classifier assigns a point the class of its single
    nearest prototype. The class of a prototype is determined by the
    position of the highest value in its soft label vector.
    
    Since we only have M prototypes, we can generate at most M distinct
    classification outcomes. However, the task requires us to correctly
    classify C different centroids, meaning we must be able to produce
    C distinct outcomes (Class 1, Class 2, ..., Class {C}).

    By the pigeonhole principle, with C required outcomes (pigeons) and at
    most M < C available outcomes (pigeonholes), at least one required
    class label cannot be produced. It is therefore impossible to
    guarantee the classification of all C centroids.
    """
    print(textwrap.dedent(explanation_part2))
    print("--> Conclusion for Part 2: Fewer than C prototypes are insufficient.\n")

    # --- Final Conclusion ---
    print("=" * 60)
    print("Final Conclusion")
    print("=" * 60)
    final_explanation = """
    Combining both parts, we've shown that C prototypes are sufficient
    and any number less than C is not. Therefore, the minimum number of
    prototypes required is exactly C.
    """
    print(textwrap.dedent(final_explanation))

    print("The relationship is expressed by the final equation:")
    print("Minimum Number of Prototypes = C")
    
    print(f"\nFor the specific case where C = {C}, the equation is:")
    print(f"Minimum Number of Prototypes = {C}")


# You can change this value to see the explanation for other cases.
number_of_classes = 5
explain_minimum_prototypes(number_of_classes)