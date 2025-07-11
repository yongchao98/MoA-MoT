import math

def calculate_max_f1_score():
    """
    Calculates the theoretically maximal F1 score for cross-language knowledge graph alignment
    under the specified constraints.

    The main constraint is that p% of relations are semantically equivalent but lexically different,
    and no external linguistic resources are available. This makes these relations un-alignable.
    """
    # p is the percentage of semantically equivalent but lexically different relations.
    p = 35.0

    # In a theoretically maximal scenario, we assume the model makes no incorrect predictions
    # for the relations it can align. Thus, False Positives (FP) = 0, and Precision = 1.0.
    precision = 1.0

    # The model cannot align the p% of lexically different relations. These become False Negatives (FN).
    # The best it can do is align the remaining (100-p)%. This percentage represents the
    # maximum possible recall.
    recall = (100.0 - p) / 100.0

    # The F1 score is the harmonic mean of precision and recall.
    # F1 = 2 * (Precision * Recall) / (Precision + Recall)
    if (precision + recall) == 0:
        f1_score = 0.0
    else:
        f1_score = 2 * (precision * recall) / (precision + recall)

    print("Problem Parameters:")
    print(f"Percentage of un-alignable relations (p): {p}%")
    print("\nDerivation of Theoretical Maximums:")
    print(f"Maximal Precision: Assumed to be {precision} (no false positives for a perfect model).")
    print(f"Maximal Recall: Limited by un-alignable relations = 1.0 - ({p}/100) = {recall}")
    print("\nF1 Score Calculation:")
    print(f"F1 = 2 * (Precision * Recall) / (Precision + Recall)")
    print(f"F1 = 2 * ({precision} * {recall}) / ({precision} + {recall})")
    print(f"F1 = {2 * precision * recall} / {precision + recall}")
    print(f"\nTheoretically Maximal F1 Score = {f1_score}")


calculate_max_f1_score()
# The final numerical answer rounded to four decimal places.
# F1 = 1.3 / 1.65 = 0.787878...
final_answer = round(1.3 / 1.65, 4)
# print(f"\nFinal Answer: <<<__{final_answer}__>>>")