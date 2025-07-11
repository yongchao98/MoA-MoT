def calculate_max_f1_score(p):
    """
    Calculates the theoretically maximal F1 score for knowledge graph alignment.

    Args:
        p (int): The percentage of relations that are semantically equivalent
                 but lexically different.
    """
    if not (0 <= p <= 100):
        print("Error: p must be a percentage between 0 and 100.")
        return

    # 1. Calculate the fraction of un-alignable and alignable relations.
    # Without external linguistic resources, p% of relations are un-alignable.
    unalignable_fraction = p / 100.0
    alignable_fraction = 1.0 - unalignable_fraction

    # 2. To achieve the maximal F1 score, we assume a perfect model for the alignable part.
    # Maximal Precision: A perfect model makes no false positive predictions.
    # Precision = TP / (TP + FP). With FP = 0, Precision = 1.
    max_precision = 1.0

    # Maximal Recall: The model can, at best, identify all the alignable relations.
    # The un-alignable relations become false negatives.
    # Recall = TP / (TP + FN) = (alignable) / (alignable + un-alignable)
    max_recall = alignable_fraction

    # 3. Calculate the F1 score, which is the harmonic mean of Precision and Recall.
    if (max_precision + max_recall) == 0:
        max_f1 = 0.0
    else:
        max_f1 = 2 * (max_precision * max_recall) / (max_precision + max_recall)

    # 4. Print the step-by-step calculation.
    print(f"Given p = {p}%, the fraction of un-alignable relations is {unalignable_fraction:.2f}.")
    print(f"Therefore, the maximum fraction of alignable relations is {alignable_fraction:.2f}.")
    print("\nTo calculate the theoretical maximum F1 score, we assume:")
    print(f" - Maximal Precision = {max_precision:.1f} (no incorrect alignments are made)")
    print(f" - Maximal Recall = {max_recall:.2f} (all alignable relations are found)")
    print("\nThe F1 score is calculated as: 2 * (Precision * Recall) / (Precision + Recall)")
    print("\nFinal Calculation:")
    print(f"Max F1 Score = 2 * ({max_precision:.1f} * {max_recall:.2f}) / ({max_precision:.1f} + {max_recall:.2f})")
    print(f"Max F1 Score = {2 * max_precision * max_recall:.2f} / {max_precision + max_recall:.2f}")
    print(f"Max F1 Score = {max_f1}")


# Given value for p
p = 35
calculate_max_f1_score(p)