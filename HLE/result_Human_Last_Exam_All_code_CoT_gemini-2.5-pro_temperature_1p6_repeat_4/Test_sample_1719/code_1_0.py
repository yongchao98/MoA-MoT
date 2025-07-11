def print_obstruction_groups(n, rank):
    """
    Prints the list of groups that classify the obstructions.

    Args:
        n: The dimension of the suspension space Sigma X, where X is a homology (n-1)-sphere.
        rank: The rank of the vector bundle, must be an even integer (2k).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(rank, int) or rank <= 0 or rank % 2 != 0:
        print("Error: rank must be a positive even integer (2k).")
        return

    print("The homotopy-theoretic obstructions for the two turnings to be homotopic are elements of a group")
    print("which is constructed from the following list of groups:\n")

    print("1. Homology groups of X:")
    print(f"   - H_{n-1}(X, Z) (The integral homology group of X in dimension n-1)")
    print(f"   - H_{n-2}(X, Z) (The integral homology group of X in dimension n-2)")
    print("")

    print("2. Homotopy groups of the structure group SO({rank}):".format(rank=rank))
    print("   - pi_1(SO({rank})) (Fundamental group)".format(rank=rank))
    print("   - pi_2(SO({rank}))".format(rank=rank))
    print("   - pi_{n}(SO({rank}))".format(n=n, rank=rank))
    print("   - pi_{n+1}(SO({rank}))".format(n=n, rank=rank))


if __name__ == '__main__':
    # You can change these values to see the groups for different dimensions.
    # n is the dimension of the base space Sigma X
    n = 4  # e.g., X is a homology 3-sphere
    # rank_2k is the rank of the vector bundle E
    rank_2k = 10 # e.g., SO(10)

    print_obstruction_groups(n, rank_2k)
    # The final answer is the conceptual list, not a specific value.
    # The list is provided by the output of the function above.
    answer_list = [
        f"H_{n-1}(X, Z)",
        f"H_{n-2}(X, Z)",
        f"pi_1(SO({rank_2k}))",
        f"pi_2(SO({rank_2k}))",
        f"pi_{n}(SO({rank_2k}))",
        f"pi_{n+1}(SO({rank_2k}))",
    ]
    # We are not providing a numerical answer, so the final format is just a list of strings
    # which represent the names of the relevant mathematical groups.
    final_answer = str(answer_list)
    # The user instruction seems to prefer a single value like <<<C>>> or <<<9.8>>>.
    # Since the answer is a list of groups, I will present it in a readable format.
    # Let's provide a compact representation for the final answer.
    final_answer_string = ', '.join(answer_list)
    print(f"\n\n<<<The obstruction groups are determined by: {final_answer_string}>>>")