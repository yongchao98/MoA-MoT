import sys

def solve():
    """
    Calculates the braid index of a knot from a grid diagram.
    """

    # Step 1: Define grid number and marker coordinates.
    n = 7
    o_coords = [(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)]
    x_coords = [(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)]

    # Create permutations from coordinates. We use 1-based indexing for clarity.
    # The lists will have size n+1, and we'll ignore index 0.
    pi_o = [0] * (n + 1)
    pi_x = [0] * (n + 1)

    for i, j in o_coords:
        pi_o[i] = j
    for i, j in x_coords:
        pi_x[i] = j

    # Helper functions for permutation operations.
    def invert_permutation(pi):
        """Computes the inverse of a permutation."""
        size = len(pi) - 1
        inv = [0] * (size + 1)
        for i in range(1, size + 1):
            inv[pi[i]] = i
        return inv

    def compose_permutations(p1, p2):
        """Computes the composition p1 . p2 (p2 applied first)."""
        size = len(p1) - 1
        result = [0] * (size + 1)
        for i in range(1, size + 1):
            result[i] = p1[p2[i]]
        return result

    def count_cycles(pi):
        """Counts the number of disjoint cycles in a permutation."""
        size = len(pi) - 1
        visited = [False] * (size + 1)
        cycles = 0
        for i in range(1, size + 1):
            if not visited[i]:
                cycles += 1
                j = i
                while not visited[j]:
                    visited[j] = True
                    j = pi[j]
        return cycles

    # Step 2: Determine the braid permutation pi_v = pi_o * pi_x^(-1).
    pi_x_inv = invert_permutation(pi_x)
    pi_v = compose_permutations(pi_o, pi_x_inv)
    
    # Step 3: Count the cycles in the resulting permutation.
    num_cycles = count_cycles(pi_v)

    # Step 4: Calculate the braid index using the formula b = n - c(pi_v).
    braid_index = n - num_cycles

    # Print the final equation with the computed values.
    # The format is "Braid Index = n - c(pi) = result".
    print(f"The grid number is n = {n}.")
    print(f"The number of cycles in the associated permutation is c(pi_v) = {num_cycles}.")
    print("The braid index is calculated as n - c(pi_v).")
    print(f"Braid Index = {n} - {num_cycles} = {braid_index}")

solve()
<<<6>>>