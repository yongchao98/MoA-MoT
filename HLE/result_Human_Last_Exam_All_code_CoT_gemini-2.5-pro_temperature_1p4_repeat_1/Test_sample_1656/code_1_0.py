import collections

def count_cycles(pi, n):
    """Counts the number of cycles in a permutation."""
    visited = [False] * (n + 1)
    cycles = 0
    for i in range(1, n + 1):
        if not visited[i]:
            cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = pi[j]
    return cycles

def solve_braid_index():
    """
    Calculates the braid index from a grid diagram's marker positions.
    """
    n = 7
    o_coords = [(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)]
    x_coords = [(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)]

    # 1. Create the pi_O and pi_X permutations
    pi_O = {col: row for col, row in o_coords}
    pi_X = {col: row for col, row in x_coords}

    # 2. Count cycles for pi_O and pi_X
    c_pi_O = count_cycles(pi_O, n)
    c_pi_X = count_cycles(pi_X, n)

    # 3. Calculate pi_H = pi_X^-1 * pi_O
    # First, find the inverse of pi_X
    pi_X_inv = {row: col for col, row in pi_X.items()}

    # Then, find the composite permutation pi_H
    pi_H = {}
    for i in range(1, n + 1):
        pi_H[i] = pi_X_inv[pi_O[i]]

    # 4. Count cycles for pi_H
    c_pi_H = count_cycles(pi_H, n)

    # 5. Calculate the braid index using the formula
    braid_index = c_pi_O + c_pi_X - c_pi_H + 1

    print("Calculation Steps:")
    print(f"Grid number n = {n}")
    print(f"Number of cycles in pi_O, c(pi_O) = {c_pi_O}")
    print(f"Number of cycles in pi_X, c(pi_X) = {c_pi_X}")
    print(f"Number of cycles in pi_H, c(pi_H) = {c_pi_H}")
    print("\nBraid Index Formula: b(K) = c(pi_O) + c(pi_X) - c(pi_H) + 1")
    print(f"Final Equation: {braid_index} = {c_pi_O} + {c_pi_X} - {c_pi_H} + 1")
    print(f"\nThe braid index of the corresponding knot is: {braid_index}")

solve_braid_index()