import sys

def solve_f1_score():
    """
    Calculates the theoretical maximal F1 score for cross-language knowledge graph alignment.

    The model assumes:
    - p: The proportion of relations that are semantically equivalent.
    - d: The Jensen-Shannon divergence, treated as a noise factor in the structural signal.
    """
    # Given parameters
    p_percent = 35
    d = 0.4

    # Convert percentage to a proportion
    p = p_percent / 100.0

    # --- Theoretical Model ---
    # The proportion of correct structural signals is (1-d).
    # The proportion of noisy/incorrect structural signals is d.

    # --- Calculating Proportions for Confusion Matrix ---
    # TP: A true equivalent exists (p) AND the signal is correct (1-d)
    tp_rate = p * (1 - d)
    # FP: No true equivalent exists (1-p) AND the signal is noisy (d)
    fp_rate = (1 - p) * d
    # FN: A true equivalent exists (p) AND the signal is noisy (d)
    fn_rate = p * d

    # --- Calculating Recall ---
    # Recall = TP / (TP + FN) = (p * (1-d)) / (p * (1-d) + p * d) = 1-d
    recall = tp_rate / (tp_rate + fn_rate)

    # --- Calculating Precision ---
    # Precision = TP / (TP + FP) = (p * (1-d)) / (p * (1-d) + (1-p) * d)
    precision = tp_rate / (tp_rate + fp_rate)

    # --- Calculating F1 Score ---
    # F1 = 2 * (Precision * Recall) / (Precision + Recall)
    f1_score = 2 * (precision * recall) / (precision + recall)

    # --- Printing the Step-by-Step Solution ---
    print("### Derivation of the Theoretical Maximal F1 Score ###")
    print("\n1. Initial Parameters:")
    print(f"  - Percentage of equivalent relations (p%): {p_percent}%")
    print(f"  - Relational distribution divergence (d): {d}")

    print("\n2. Calculating Theoretical Recall:")
    print("  - Formula: Recall = (p * (1-d)) / (p * (1-d) + p * d) = 1 - d")
    print(f"  - Calculation: Recall = 1 - {d} = {recall:.4f}")

    print("\n3. Calculating Theoretical Precision:")
    print("  - Formula: Precision = (p * (1-d)) / (p * (1-d) + (1-p) * d)")
    print(f"  - Numerator (TP rate): {p:.2f} * (1 - {d:.2f}) = {tp_rate:.4f}")
    print(f"  - Denominator (TP+FP rate): ({p:.2f} * (1 - {d:.2f})) + ((1 - {p:.2f}) * {d:.2f}) = {tp_rate + fp_rate:.4f}")
    print(f"  - Calculation: Precision = {tp_rate:.4f} / {tp_rate + fp_rate:.4f} = {precision:.4f}")

    print("\n4. Calculating the Final F1 Score:")
    print("  - Formula: F1 = 2 * (Precision * Recall) / (Precision + Recall)")
    print(f"  - Calculation: F1 = 2 * ({precision:.4f} * {recall:.4f}) / ({precision:.4f} + {recall:.4f})")
    print(f"  - Calculation: F1 = {2 * precision * recall:.4f} / {precision + recall:.4f}")
    
    print("\n" + "="*50)
    print(f"The theoretically maximal F1 score is: {f1_score:.4f}")
    print("="*50)

    # Suppress final answer line in stdout to not confuse the user
    original_stdout = sys.stdout
    sys.stdout = open('/dev/null', 'w')
    print(f"<<<{f1_score:.3f}>>>")
    sys.stdout = original_stdout

solve_f1_score()