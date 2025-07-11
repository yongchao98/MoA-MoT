import sys

def solve():
    """
    Calculates the theoretically maximal F1 score for cross-language knowledge graph
    alignment given a percentage of lexically different relations and no external
    linguistic resources.
    """
    # p: The percentage of semantically equivalent but lexically different relations.
    p = 35

    # In a theoretical best-case scenario, the model makes no incorrect predictions (False Positives = 0).
    # Precision = True Positives / (True Positives + False Positives)
    # With FP = 0, Precision becomes 1.0.
    precision = 1.0

    # Without external linguistic resources, the model cannot identify the p% of
    # relations that are lexically different. These are missed alignments (False Negatives).
    # The model can only find the (100-p)% of relations that are lexically accessible.
    # This directly determines the maximum recall.
    # Recall = True Positives / (True Positives + False Negatives) = (100 - p) / 100
    recall = (100.0 - p) / 100.0

    # The F1 score is the harmonic mean of precision and recall.
    # F1 = 2 * (Precision * Recall) / (Precision + Recall)
    if (precision + recall) == 0:
        f1_score = 0.0
    else:
        f1_score = 2 * (precision * recall) / (precision + recall)

    print(f"Calculating the theoretical maximum F1 score with p = {p}%...\n")
    print("Step 1: Determine Precision.")
    print("A perfect model makes no mistakes, so False Positives = 0.")
    print(f"Precision = TP / (TP + 0) = 1.0\n")

    print("Step 2: Determine Recall.")
    print("Without linguistic resources, the model cannot align the 35% of lexically different relations.")
    print("Recall = (100 - p) / 100")
    print(f"Recall = (100 - {p}) / 100 = {recall}\n")

    print("Step 3: Calculate F1 Score.")
    print("F1 Score = 2 * (Precision * Recall) / (Precision + Recall)")
    # Final Equation Output
    print(f"F1 Score = 2 * ({precision} * {recall}) / ({precision} + {recall})")
    print(f"F1 Score = {2 * precision * recall} / {precision + recall}")
    print(f"\nThe theoretically maximal F1 score is: {f1_score}")
    
    # Required for automated checking
    sys.stdout.write(f"<<<{f1_score}>>>")

solve()