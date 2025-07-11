def calculate_mistake_bound(n, c):
    """
    Calculates an upper bound on the number of mistakes made by a majority
    voting algorithm in a specific variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert must make to be removed.
                 The true expert makes strictly fewer than c mistakes.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of experts 'n' must be a positive integer.")
        return
    if not isinstance(c, int) or c <= 0:
        print("Error: The mistake cutoff 'c' must be a positive integer.")
        return
    if n == 1 and c == 1:
        # A single true expert can make 0 mistakes, so the algorithm makes 0 mistakes.
        bound = 0
        k_A_bound = 0
        k_B_bound = 0
        print(f"For n = {n} and c = {c}:")
        print("With only one expert who is always right (makes < 1 mistakes), the algorithm also makes 0 mistakes.")
        print(f"Bound = {bound}")
        return


    # The bound on Type A mistakes (k_A) is (n-1) * c
    k_A_bound = (n - 1) * c
    # The bound on Type B mistakes (k_B) is c - 1
    k_B_bound = c - 1

    # The total bound is the sum of the two bounds
    total_bound = k_A_bound + k_B_bound

    print(f"An upper bound on the number of algorithm mistakes can be calculated as follows:")
    print(f"Let n = {n} (number of experts) and c = {c} (mistake cutoff).")
    print("\nThe total mistakes (M) can be split into two types:")
    print("1. (k_A): Algorithm is wrong, true expert is right.")
    print("2. (k_B): Algorithm is wrong, true expert is wrong.")
    print("\nBounding k_A:")
    print(f"The (n-1) = {n-1} false experts can make at most c = {c} mistakes each.")
    print(f"Total mistake capacity for false experts = ({n} - 1) * {c} = {n-1} * {c} = {k_A_bound}.")
    print(f"Each Type A mistake requires at least one false expert to be wrong, so k_A <= {k_A_bound}.")
    print("\nBounding k_B:")
    print(f"The true expert makes strictly fewer than c mistakes, so at most c-1 mistakes.")
    print(f"So, k_B <= {c} - 1 = {k_B_bound}.")
    print("\nTotal Bound Calculation:")
    print(f"M = k_A + k_B <= ({n} - 1) * {c} + ({c} - 1)")
    print(f"M <= {n-1} * {c} + {k_B_bound}")
    print(f"M <= {k_A_bound} + {k_B_bound}")
    print(f"M <= {total_bound}")
    print(f"\nFinal simplified upper bound is n*c - 1 = {n}*{c} - 1 = {n*c-1}")


# --- Example Usage ---
# You can change these values to see the bound for different scenarios.
example_n = 10
example_c = 5

calculate_mistake_bound(example_n, example_c)