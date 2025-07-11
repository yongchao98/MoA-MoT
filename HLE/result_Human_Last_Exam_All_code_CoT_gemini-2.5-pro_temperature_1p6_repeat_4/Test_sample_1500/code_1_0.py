import math

def calculate_b2_sun_orbit(n, partition):
    """
    Calculates the second Betti number for a coadjoint orbit of SU(n).

    Args:
        n (int): The rank of the group SU(n).
        partition (list[int]): A list of integers representing the sizes of blocks
                               of equal eigenvalues for lambda. The sum of the
                               partition must be equal to n.
                               For example, for SU(4):
                               - Regular orbit: [1, 1, 1, 1]
                               - Singular orbit (like CP^2 grassmannian): [2, 2]
                               - Trivial orbit: [4]

    Returns:
        tuple: A tuple containing (b_2, n-1, is_match).
    """
    if sum(partition) != n:
        raise ValueError(f"The sum of the partition elements {partition} must be equal to n={n}")

    # The number of distinct eigenvalues k is the length of the partition.
    k = len(partition)

    # For a coadjoint orbit of SU(n) corresponding to a partition of n,
    # the second Betti number b_2 is k-1, where k is the number of parts.
    b2_actual = k - 1
    
    b2_conjecture = n - 1
    
    return b2_actual, b2_conjecture, b2_actual == b2_conjecture

def demonstrate_question_b():
    """
    Demonstrates that b_2 is not always n-1 for orbits in SU(n).
    """
    print("--- Demonstrating Question (b) for G = SU(n) ---")
    print("Is the second Betti number b_2(O_lambda) always given by n-1?\n")

    # Case 1: SU(3)
    n = 3
    print(f"--- Group SU({n}) ---")
    print(f"The conjectured value for b_2 is n-1 = {n-1}\n")

    # Regular orbit (lambda has distinct eigenvalues)
    partition_regular = [1, 1, 1]
    b2, conj, match = calculate_b2_sun_orbit(n, partition_regular)
    print(f"For a REGULAR orbit (partition {partition_regular}):")
    print(f"The number of distinct eigenvalues is k = {len(partition_regular)}.")
    print(f"The actual b_2 is k-1 = {b2}.")
    print(f"Does it match the conjectured value? {match}\n")

    # Singular orbit (e.g., leads to CP^2)
    partition_singular = [2, 1]
    b2, conj, match = calculate_b2_sun_orbit(n, partition_singular)
    print(f"For a SINGULAR orbit (partition {partition_singular}):")
    print(f"The number of distinct eigenvalues is k = {len(partition_singular)}.")
    print(f"The actual b_2 is k-1 = {b2}.")
    print(f"Does it match the conjectured value? {match}\n")

    # Trivial orbit (lambda = 0)
    partition_trivial = [3]
    b2, conj, match = calculate_b2_sun_orbit(n, partition_trivial)
    print(f"For the TRIVIAL orbit (partition {partition_trivial}):")
    print(f"The number of distinct eigenvalues is k = {len(partition_trivial)}.")
    print(f"The actual b_2 is k-1 = {b2}.")
    print(f"Does it match the conjectured value? {match}\n")
    
    print("Conclusion: Since singular orbits provide counterexamples, the statement is FALSE.")


if __name__ == '__main__':
    demonstrate_question_b()