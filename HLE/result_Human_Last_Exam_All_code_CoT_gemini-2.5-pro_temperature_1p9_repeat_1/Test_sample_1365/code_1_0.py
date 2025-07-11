import math

def calculate_mistake_bound(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes made
    by a majority voting algorithm under specific expert problem conditions.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes after which a false expert is removed.
                 The true expert makes strictly fewer than c mistakes.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be positive integers.")
        return

    # Bound for M_w: mistakes when the true expert is wrong.
    # M_w <= K_true <= c - 1
    m_w_bound = c - 1

    # Bound for M_c: mistakes when the true expert is correct.
    # M_c <= (n - 1) * c / 2
    m_c_bound = (n - 1) * c / 2

    # Total bound M = M_w + M_c
    total_bound = m_w_bound + m_c_bound
    
    # The number of mistakes M must be an integer.
    final_integer_bound = math.floor(total_bound)

    print("Derivation of the upper bound for M (total mistakes):")
    print(f"Number of experts (n) = {n}")
    print(f"Mistake tolerance (c) = {c}\n")

    print("1. The total mistakes M can be split into M_w and M_c.")
    print("   M = M_w (true expert is wrong) + M_c (true expert is correct).\n")

    print("2. Bounding M_w:")
    print(f"   M_w <= K_true (mistakes by true expert) <= c - 1")
    print(f"   M_w <= {c} - 1 = {m_w_bound}\n")
    
    print("3. Bounding M_c:")
    print(f"   M_c <= (n - 1) * c / 2")
    print(f"   M_c <= ({n} - 1) * {c} / 2 = {n-1} * {c} / 2 = {m_c_bound}\n")

    print("4. Combining the bounds for the total mistakes M:")
    print(f"   M <= M_w_bound + M_c_bound")
    print(f"   M <= {m_w_bound} + {m_c_bound} = {total_bound}\n")
    
    print("Since the number of mistakes must be an integer:")
    print(f"Final Answer: M <= {final_integer_bound}")


if __name__ == '__main__':
    try:
        # Example usage with user input
        n_input = int(input("Enter the number of experts (n): "))
        c_input = int(input("Enter the mistake tolerance (c): "))
        print("-" * 30)
        calculate_mistake_bound(n_input, c_input)
    except ValueError:
        print("Invalid input. Please enter integers for n and c.")
