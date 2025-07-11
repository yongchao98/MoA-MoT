import sys

def solve():
    """
    Calculates the theoretically maximal F1 score for a cross-language KG alignment task
    under the specified conditions.
    """
    # p is the percentage of semantically equivalent but lexically different relations.
    p = 35.0

    print(f"Analyzing the theoretical maximum F1 score with p = {p}%.\n")

    # In a theoretically maximal scenario, an embedding-based model without external
    # linguistic resources cannot bridge the lexical gap for the p% of relations.
    # Therefore, these relations are considered un-alignable.
    # The best possible model can perfectly align the remaining (100-p)% of relations.

    # We normalize the total number of true alignments to 1 for simplicity.
    total_positive_relations = 1.0

    # True Positives (TP): The proportion of relations that are correctly aligned.
    # This corresponds to the relations that are alignable.
    tp = total_positive_relations * (100.0 - p) / 100.0

    # False Negatives (FN): The proportion of relations that are missed.
    # This corresponds to the un-alignable relations.
    fn = total_positive_relations * p / 100.0

    # False Positives (FP): A "theoretically maximal" model is assumed to be
    # perfectly precise, so it makes no incorrect alignments.
    fp = 0.0

    # Step 1: Calculate Precision
    # Precision = TP / (TP + FP)
    # This measures the accuracy of the alignments made.
    if (tp + fp) == 0:
      precision = 0.0
    else:
      precision = tp / (tp + fp)

    # Step 2: Calculate Recall
    # Recall = TP / (TP + FN)
    # This measures how many of the truly positive alignments were found.
    if (tp + fn) == 0:
        recall = 0.0
    else:
        recall = tp / (tp + fn)

    # Step 3: Calculate the F1 Score
    # F1 Score is the harmonic mean of Precision and Recall.
    if (precision + recall) == 0:
        f1_score = 0.0
    else:
        f1_score = 2 * (precision * recall) / (precision + recall)

    # Print the detailed breakdown of the calculation as requested.
    print("Calculation Breakdown:")
    print(f"Precision = TP / (TP + FP) = {tp:.2f} / ({tp:.2f} + {fp:.2f}) = {precision:.4f}")
    print(f"Recall    = TP / (TP + FN) = {tp:.2f} / ({tp:.2f} + {fn:.2f}) = {recall:.4f}")
    print("\nThe final F1 score is calculated using the formula:")
    print("F1 = 2 * (Precision * Recall) / (Precision + Recall)")
    print(f"F1 = 2 * ({precision:.4f} * {recall:.4f}) / ({precision:.4f} + {recall:.4f})")
    print(f"F1 = {2 * precision * recall:.4f} / {precision + recall:.4f}")
    print(f"Maximal F1 Score = {f1_score:.4f}")

solve()