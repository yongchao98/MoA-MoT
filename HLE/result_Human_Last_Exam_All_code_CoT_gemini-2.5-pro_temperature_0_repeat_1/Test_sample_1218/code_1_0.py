def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the derived formula.

    The problem asks for the maximum value of n such that there exists a
    k-uniform intersecting family F with full differences of size k-1.
    
    1. The condition that every (k-1)-subset S is a difference F \\ F' implies
       that for any S, there exist F, F' in the family such that |F intersect F'| = 1.
    2. Let F = S U {x} and F' = {x} U T. For F' to exist, the set from which T
       is drawn, [n] \\ F, must have at least k-1 elements.
    3. This leads to the inequality n - k >= k - 1, which simplifies to n >= 2k - 1.
    4. We can show that n = 2k - 1 is achievable by taking the family F to be
       all k-subsets of a (2k-1)-element set.
    5. It can be proven that for n >= 2k, no such family exists.
    
    Thus, the maximum value of n is 2k - 1.
    """
    if not isinstance(k, int) or k < 2:
        print("k must be an integer greater than or equal to 2.")
        return

    # Calculate n using the formula n = 2k - 1
    n = 2 * k - 1
    
    # The final equation is n = 2 * k - 1
    print(f"For k = {k}, the maximum value of n is determined by the equation:")
    print(f"n = 2 * {k} - 1 = {n}")
    
    # As requested, output each number in the final equation
    print("\nThe numbers in the final equation are:")
    print(2)
    print(k)
    print(1)
    print(n)

# Example usage with a specific value for k.
# For the purpose of this demonstration, let's use k=4.
k_example = 4
solve_max_n(k_example)