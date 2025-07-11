import math

# --- Configuration ---
# You can change these values to see the bound for different scenarios.
# n: the total number of experts.
n = 10
# c: the number of mistakes an expert must make to be removed.
c = 5
# -------------------

if not isinstance(n, int) or not isinstance(c, int) or n <= 1 or c < 1:
    print("Error: Please provide integer values for n and c, where n > 1 and c >= 1.")
else:
    # Part 1: Bound on mistakes when the true expert is wrong.
    # The true expert makes at most c-1 mistakes.
    m_w_bound = c - 1

    # Part 2: Bound on mistakes when the true expert is right.
    # This is bounded by the total number of mistakes the other n-1 experts can make.
    # Each of the n-1 experts can make at most c mistakes.
    m_r_bound = (n - 1) * c

    # The total upper bound is the sum of these two parts.
    upper_bound = m_w_bound + m_r_bound

    # --- Output ---
    print(f"Given n = {n} experts and a mistake threshold c = {c}:")
    print("\n1. The number of mistakes when the true expert is wrong is at most c - 1.")
    print(f"   Calculation: {c} - 1 = {m_w_bound}")

    print("\n2. The number of mistakes when the true expert is right is at most (n - 1) * c.")
    print(f"   Calculation: ({n} - 1) * {c} = {m_r_bound}")

    print("\n3. The total upper bound on mistakes is the sum of these two values.")
    print("   Total Mistakes <= (Bound from #1) + (Bound from #2)")
    print(f"\nFinal Equation:")
    print(f"({c} - 1) + ({n} - 1) * {c} = {upper_bound}")
    # You can also simplify the formula to n*c - 1
    # print(f"Simplified: {n} * {c} - 1 = {upper_bound}")

    # The final answer required by the problem format
    final_answer_value = upper_bound
    # print(f"\n<<< {final_answer_value} >>>")