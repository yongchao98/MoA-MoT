def identify_obstructions(n, rank):
    """
    Identifies and prints the homotopy-theoretic obstructions for the given problem.

    Args:
        n (int): An integer determined by the homology of X, where H_{n-1}(X) is the
                 only non-trivial positive-dimensional homology group.
        rank (int): The rank of the vector bundle E, which must be an even integer (2k).
    """
    if rank % 2 != 0:
        print("Error: The rank of the vector bundle must be an even integer (2k).")
        return

    k = rank // 2

    print("The homotopy-theoretic obstructions for the two paths to be homotopic are determined by an element")
    print("in the group pi_1(Aut(E)). This group is an extension of a quotient of one homotopy group by a")
    print("subgroup of another. The groups that classify these obstructions are listed below.")
    print("-" * 30)
    
    # The list of groups that determine the obstruction.
    # The first one defines the integer n.
    print(f"1. H_{n-1}(X): The homology group of X that defines the integer n. By definition, H_{n-1}(X) = Z.")
    
    # The other groups are homotopy groups of SO(2k).
    print(f"2. pi_1(SO({rank})): The fundamental group of the special orthogonal group SO({rank}).")
    print(f"3. pi_{n}(SO({rank})): The {n}-th homotopy group of SO({rank}). The index {n} is derived from the dimension of X.")
    print(f"4. pi_{n+1}(SO({rank})): The {n+1}-th homotopy group of SO({rank}). The index {n+1} is also derived from the dimension of X.")
    print("-" * 30)


if __name__ == '__main__':
    # --- Example Case ---
    # Let X be a homology 3-sphere. This means n-1=3, so n=4.
    example_n = 4
    # Let the vector bundle E be of rank 10. So 2k=10.
    example_rank = 10
    
    print(f"Running for the example case: n = {example_n}, rank = {example_rank}")
    identify_obstructions(example_n, example_rank)
