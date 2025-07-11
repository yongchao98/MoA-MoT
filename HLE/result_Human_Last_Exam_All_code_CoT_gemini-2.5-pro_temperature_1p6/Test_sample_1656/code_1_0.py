def solve():
    """
    Calculates the braid index of a knot from a grid diagram.
    """
    n = 7
    # O_pos[column] = row
    O_pos = {1:1, 2:7, 3:4, 4:5, 5:3, 6:6, 7:2}
    # X_pos[column] = row
    X_pos = {1:2, 2:6, 3:3, 4:1, 5:7, 6:5, 7:4}

    # Helper maps for easier permutation calculation
    # O_row_to_col[row] = column
    O_row_to_col = {v: k for k, v in O_pos.items()}
    # X_row_to_col[row] = column
    X_row_to_col = {v: k for k, v in X_pos.items()}

    # 1. Calculate the column permutation (pi_col)
    # pi_col(i) = j: O in col i is at row k, X in row k is at col j.
    # k = O_pos[i] -> j = X_row_to_col[k]
    pi_col = {}
    for i in range(1, n + 1):
        k = O_pos[i]
        pi_col[i] = X_row_to_col[k]

    # 2. Calculate the row permutation (pi_row)
    # pi_row(i) = j: O in row i is at col k, X in col k is at row j.
    # k = O_row_to_col[i] -> j = X_pos[k]
    pi_row = {}
    for i in range(1, n + 1):
        k = O_row_to_col[i]
        pi_row[i] = X_pos[k]

    def count_cycles(perm, size):
        """Counts the number of cycles in a permutation."""
        visited = [False] * (size + 1)
        cycles = 0
        for i in range(1, size + 1):
            if not visited[i]:
                cycles += 1
                j = i
                while not visited[j]:
                    visited[j] = True
                    j = perm[j]
        return cycles

    # 3. Count cycles in each permutation
    n_col = n
    c_pi_col = count_cycles(pi_col, n)
    c_pi_row = count_cycles(pi_row, n)

    # 4. Apply the formula for the braid index
    braid_index = n_col - c_pi_col + c_pi_row

    # Output the final calculation step-by-step
    print("The braid index is calculated using the formula:")
    print("Braid Index = n - c(pi_col) + c(pi_row)")
    print(f"where n is the grid number, c(pi_col) is the number of cycles in the column permutation, and c(pi_row) is the number of cycles in the row permutation.")
    print(f"\nCalculated values:")
    print(f"n = {n_col}")
    print(f"c(pi_col) = {c_pi_col}")
    print(f"c(pi_row) = {c_pi_row}")
    print(f"\nFinal calculation:")
    print(f"Braid Index = {n_col} - {c_pi_col} + {c_pi_row} = {braid_index}")
    print(f"\nThe braid index of the corresponding knot is {braid_index}.")


solve()
<<<7>>>