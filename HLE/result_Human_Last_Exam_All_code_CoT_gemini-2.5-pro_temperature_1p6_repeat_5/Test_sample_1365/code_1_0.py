import math

def calculate_mistake_bound(n, c):
    """
    Calculates and explains the upper bound on mistakes for a variant of the experts problem.

    Args:
        n (int): The number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Please provide positive integers for n and c.")
        return

    # Step 1: Explain the setup
    print("Derivation of the upper bound on algorithm mistakes (M):")
    print("-" * 50)
    print(f"We have n = {n} experts and a mistake removal threshold c = {c}.")
    print("One 'true expert' makes strictly fewer than c mistakes.")
    print("Other experts are removed after making c mistakes.")
    print("-" * 50)

    # Step 2: Bound the total mistakes made by all experts combined
    total_expert_mistakes_bound = (n - 1) * c + (c - 1)
    print("1. Find an upper bound for the total mistakes made by all experts (C_total).")
    print("   - The 'true expert' makes at most (c - 1) mistakes.")
    print(f"     c - 1 = {c} - 1 = {c - 1}")
    print("   - The other (n - 1) experts make at most c mistakes each.")
    print(f"     (n - 1) * c = ({n} - 1) * {c} = {(n - 1) * c}")
    print("   - Therefore, the total number of expert mistakes is bounded:")
    print(f"     C_total <= (n - 1) * c + (c - 1)")
    print(f"     C_total <= ({n} - 1) * {c} + ({c} - 1) = {total_expert_mistakes_bound}")
    print("-" * 50)

    # Step 3: Relate algorithm mistakes (M) to total expert mistakes (C_total)
    print("2. Relate algorithm mistakes (M) to C_total.")
    print("   - When the algorithm makes a mistake, more than half of the active experts must have also been wrong.")
    print("   - The true expert is never removed, so the number of active experts is always at least 1.")
    print("   - Let's sum this effect over all M mistakes the algorithm makes.")
    print("   - Each mistake adds more than 1/2 * (number of active experts) to C_total.")
    print("   - At a minimum, C_total increases by more than 1/2 for each algorithm mistake.")
    print("   - This gives the inequality: C_total > M / 2.")
    print("-" * 50)

    # Step 4: Combine the bounds to find the upper bound for M
    strict_bound = 2 * (n * c - 1)
    integer_bound = strict_bound - 1
    print("3. Combine the two inequalities to find the bound on M.")
    print(f"   From step 1: C_total <= {total_expert_mistakes_bound}")
    print(f"   From step 2: C_total > M / 2")
    print(f"   Combining them: {total_expert_mistakes_bound} >= C_total > M / 2")
    print(f"   So: {total_expert_mistakes_bound} > M / 2")
    print(f"   Multiplying by 2: 2 * {total_expert_mistakes_bound} > M")
    print(f"   Which means: {strict_bound} > M")
    print("-" * 50)
    print("Final Result:")
    print(f"The number of mistakes M must be an integer strictly less than {strict_bound}.")
    print(f"So, a valid upper bound for M is {strict_bound} - 1 = {integer_bound}.")

# Example usage with n=10 experts and c=5 as the mistake threshold
calculate_mistake_bound(n=10, c=5)