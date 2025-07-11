def calculate_max_f1_score(p_percent):
    """
    Calculates the theoretically maximal F1 score given the percentage of
    semantically equivalent but lexically different relations.

    Args:
        p_percent (float): The percentage of relations that are lexically different.
    """
    if not 0 <= p_percent <= 100:
        print("Error: p must be between 0 and 100.")
        return

    p = p_percent / 100.0

    # In a "theoretically maximal" scenario, the model makes no incorrect predictions.
    # TP = All alignable relations, FP = 0.
    # Precision = TP / (TP + FP) = TP / TP = 1
    precision = 1.0

    # Recall is limited by the relations the model cannot possibly identify.
    # The model can only identify (1 - p) of the total relations.
    # Recall = TP / (TP + FN) = (1 - p) / ((1 - p) + p) = 1 - p
    recall = 1.0 - p

    # F1 score is the harmonic mean of precision and recall.
    if (precision + recall) == 0:
        f1_score = 0.0
    else:
        f1_score = 2 * (precision * recall) / (precision + recall)

    # For the final output, we need to show the numbers in the equation.
    numerator = 2 * precision * recall
    denominator = precision + recall

    print(f"Given p = {p_percent}% lexically different relations:")
    print("--------------------------------------------------")
    print("1. We assume the model cannot align these relations without external resources.")
    print(f"2. Therefore, the maximum possible Recall is (1 - {p}) = {recall:.2f}")
    print("3. For a 'maximal' score, we assume the model makes no incorrect alignments, so Precision is {precision:.1f}")
    print("\nCalculating the F1 Score:")
    print(f"   F1 = 2 * (Precision * Recall) / (Precision + Recall)")
    print(f"   F1 = 2 * ({precision:.1f} * {recall:.2f}) / ({precision:.1f} + {recall:.2f})")
    print(f"   F1 = {numerator:.2f} / {denominator:.2f}")
    print(f"\nTheoretically maximal F1 score: {f1_score:.4f}")

# The value of p from the problem description
p_value = 35
calculate_max_f1_score(p_value)

# Final answer derivation for direct output
p_final = p_value / 100.0
recall_final = 1.0 - p_final
precision_final = 1.0
f1_final = 2 * (precision_final * recall_final) / (precision_final + recall_final)
# print(f"<<<{f1_final:.4f}>>>")