def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by a majority voting algorithm
    in a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert must make to be removed.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of experts) must be a positive integer.")
        return
    if not isinstance(c, int) or c <= 0:
        print("Error: c (mistake threshold) must be a positive integer.")
        return

    # The upper bound is derived as nc - 1.
    upper_bound = n * c - 1

    print(f"For a problem with n = {n} experts and a mistake threshold c = {c}:")
    print("The upper bound on algorithm mistakes is given by the formula: n * c - 1")
    # We explicitly print each number in the final equation as requested.
    print(f"Bound = ({n} * {c}) - 1 = {upper_bound}")

# Example usage with some values for n and c.
num_experts = 10
mistake_threshold = 5
calculate_mistake_bound(num_experts, mistake_threshold)