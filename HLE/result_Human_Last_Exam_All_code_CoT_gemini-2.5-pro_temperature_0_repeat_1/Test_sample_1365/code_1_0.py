import argparse

def calculate_mistake_bound(n: int, c: int):
    """
    Calculates and explains the upper bound on mistakes for a variant of the experts problem.

    Args:
        n: The total number of experts.
        c: The mistake threshold for removing an expert.
    """
    if n < 1 or c < 1:
        print("Error: Number of experts (n) and mistake threshold (c) must be at least 1.")
        return

    # M1 is bounded by the number of mistakes of the true expert
    m1_bound = c - 1

    # M2 is bounded by half the total possible mistakes of the false experts
    m2_bound = (n - 1) * c / 2

    # Total bound
    total_bound = m1_bound + m2_bound

    print("An upper bound on the number of mistakes is derived as follows:")
    print("Let M be the total mistakes. M = M1 + M2, where:")
    print("  - M1: Algorithm mistakes when the true expert is wrong.")
    print("  - M2: Algorithm mistakes when the true expert is right.")
    print("\n--- Step 1: Bounding M1 ---")
    print(f"The true expert makes strictly fewer than c mistakes (m_true < {c}).")
    print(f"So, M1 <= m_true <= c - 1")
    print(f"M1 <= {c} - 1 = {m1_bound}")

    print("\n--- Step 2: Bounding M2 ---")
    print(f"The (n-1) false experts can make at most (n-1)*c mistakes in total.")
    print(f"Total false expert mistakes <= ({n} - 1) * {c} = {(n - 1) * c}")
    print("For each M2 mistake, at least 2 false experts must be wrong.")
    print(f"This means 2 * M2 <= Total false expert mistakes.")
    print(f"So, M2 <= (n-1)*c / 2")
    print(f"M2 <= ({n} - 1) * {c} / 2 = {(n-1)*c}/2 = {m2_bound}")

    print("\n--- Step 3: Final Upper Bound ---")
    print("M = M1 + M2 <= (c - 1) + (n - 1) * c / 2")
    print(f"M <= ({c} - 1) + ({n} - 1) * {c} / 2")
    print(f"M <= {m1_bound} + {m2_bound}")
    print(f"Final Upper Bound = {total_bound}")


if __name__ == '__main__':
    # Set up argument parser to get n and c from the command line
    parser = argparse.ArgumentParser(
        description="Calculate the upper bound on mistakes for a majority voting algorithm.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-n',
        type=int,
        required=True,
        help='The total number of experts.'
    )
    parser.add_argument(
        '-c',
        type=int,
        required=True,
        help='The mistake threshold for removing an expert.'
    )
    
    # Example usage:
    # python your_script_name.py -n 101 -c 10

    args = parser.parse_args()
    calculate_mistake_bound(args.n, args.c)