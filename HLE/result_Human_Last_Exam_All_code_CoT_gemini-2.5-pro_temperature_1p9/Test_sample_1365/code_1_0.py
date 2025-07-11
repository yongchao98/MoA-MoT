import sys

def solve_experts_problem(n, c):
    """
    Calculates an upper bound on the number of mistakes for a variant of the experts problem.

    Args:
      n: The number of experts.
      c: The number of mistakes after which an expert is removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 0 or c <= 0:
        print("Error: n and c must be positive integers.")
        return

    # M2 is bounded by the number of mistakes the true expert can make.
    # The true expert makes strictly fewer than c mistakes, so at most c-1.
    m2_bound = c - 1

    # M1 is bounded by the total mistakes the other n-1 experts can make.
    # Each can make at most c mistakes before being removed.
    n_minus_1 = n - 1
    m1_bound = n_minus_1 * c
    
    # The total bound is the sum of the bounds for M1 and M2.
    total_bound = m1_bound + m2_bound

    print("Deriving the upper bound:")
    print(f"Let M be the total number of mistakes by the algorithm.")
    print(f"We split M into M1 (when the true expert is right) and M2 (when the true expert is wrong).")
    print(f"The number of mistakes the true expert makes is less than c, so M2 <= c - 1.")
    print(f"For n={n} and c={c}, M2 <= {c} - 1 = {c_minus_1}.")
    print(f"The other n-1 experts can make at most (n-1)*c mistakes in total.")
    print(f"Each mistake in M1 requires at least one of these n-1 experts to be wrong.")
    print(f"Therefore, M1 <= (n-1) * c = ({n}-1) * {c} = {m1_bound}.")
    print(f"The final upper bound M = M1 + M2 is calculated as:")
    print(f"M <= ({n} - 1) * {c} + ({c} - 1) = {n_minus_1} * {c} + {c_minus_1} = {m1_bound} + {m2_bound} = {total_bound}")


# Example usage with sample values for n and c.
# You can change these values to test with other numbers.
if __name__ == '__main__':
    # You can provide n and c as command-line arguments, e.g., python your_script.py 10 5
    if len(sys.argv) == 3:
        try:
            n_experts = int(sys.argv[1])
            c_mistakes = int(sys.argv[2])
            solve_experts_problem(n_experts, c_mistakes)
        except ValueError:
            print("Please provide integer values for n and c.")
    else:
        # Default example values if no command-line arguments are given
        print("Running with default values n=10, c=5.")
        print("You can also run as 'python your_script.py <n> <c>' to provide your own values.")
        n_experts = 10
        c_mistakes = 5
        solve_experts_problem(n_experts, c_mistakes)