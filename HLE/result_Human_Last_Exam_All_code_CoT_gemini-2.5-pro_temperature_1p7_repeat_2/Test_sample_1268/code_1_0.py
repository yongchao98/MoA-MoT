import math

def calculate_minkowski_bound(N):
    """
    Calculates the Minkowski bound for the real quadratic field Q(sqrt(N)).

    Args:
        N: A squarefree natural number.
    """
    print(f"Analyzing the field Q(sqrt(N)) for N = {N}\n")

    # Step 1: Determine the discriminant of the field.
    # The discriminant (Δ) depends on whether N is congruent to 1, 2, or 3 modulo 4.
    if N % 4 == 1:
        delta = N
        case_explanation = f"Since N = {N} is congruent to 1 (mod 4), the discriminant is Δ = N."
    else:  # Covers cases where N % 4 is 2 or 3
        delta = 4 * N
        case_explanation = f"Since N = {N} is congruent to 2 or 3 (mod 4), the discriminant is Δ = 4 * N."

    print(f"Step 1: Determine the discriminant (Δ).")
    print(case_explanation)
    if N % 4 != 1:
      print(f"Δ = 4 * {N} = {delta}")
    else:
      print(f"Δ = {delta}")


    # Step 2: Calculate the Minkowski bound (M_K).
    # For a real quadratic field, the formula is M_K = (1/2) * sqrt(|Δ|).
    minkowski_bound = 0.5 * math.sqrt(abs(delta))

    print("\nStep 2: Calculate the Minkowski bound (M_K).")
    print("The formula for the upper bound is: M_K = (1/2) * sqrt(|Δ|)")
    print("Plugging in the value for Δ:")
    # Print the equation with the specific numbers, as requested.
    print(f"M_K = {0.5} * sqrt({delta})")
    print(f"The resulting upper bound is: M_K = {minkowski_bound}")

# --- Main execution ---
# Let's use a squarefree number, for example, N = 105 (3 * 5 * 7).
# 105 mod 4 = 1.
squarefree_N = 105
calculate_minkowski_bound(squarefree_N)
print("-" * 20)
# Let's use another squarefree number, N = 30 (2 * 3 * 5).
# 30 mod 4 = 2.
squarefree_N_2 = 30
calculate_minkowski_bound(squarefree_N_2)
