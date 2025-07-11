import collections

def solve_braid_index():
    """
    Calculates the braid index from a grid diagram.
    """
    # Grid size
    n = 7

    # Positions of 'o' and 'x' marks using 1-based indexing for columns and rows.
    o_positions = [(1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)]
    x_positions = [(1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)]

    # We use 0-based indexing for lists in Python.
    # pi_O[i] = j means 'o' is in column i+1 and row j+1.
    pi_O = [0] * n
    for col, row in o_positions:
        pi_O[col - 1] = row - 1

    # pi_X[i] = j means 'x' is in column i+1 and row j+1.
    pi_X = [0] * n
    for col, row in x_positions:
        pi_X[col - 1] = row - 1
        
    # To compute alpha = pi_X_inv * pi_O, we first need the inverse of pi_X.
    # If pi_X[i] = j, then pi_X_inv[j] = i.
    pi_X_inv = [0] * n
    for i in range(n):
        pi_X_inv[pi_X[i]] = i

    # Now, compose the permutations to get alpha.
    # alpha[i] = pi_X_inv(pi_O(i))
    alpha = [0] * n
    for i in range(n):
        alpha[i] = pi_X_inv[pi_O[i]]

    # Next, we count the number of cycles in the permutation alpha.
    visited = [False] * n
    num_cycles = 0
    for i in range(n):
        if not visited[i]:
            num_cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = alpha[j]

    c_alpha = num_cycles

    # The braid index is given by the formula b = n - c(alpha).
    braid_index = n - c_alpha

    print("The braid index of a knot from a grid diagram is calculated using the formula: b = n - c(alpha)")
    print("Where 'n' is the grid number and 'c(alpha)' is the number of cycles in the row-action permutation.")
    print(f"For the given grid:")
    print(f"Grid number n = {n}")
    print(f"Number of cycles c(alpha) = {c_alpha}")
    print("The braid index is calculated as:")
    print(f"b = {n} - {c_alpha} = {braid_index}")

solve_braid_index()