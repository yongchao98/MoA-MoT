def solve_scl_problem():
    """
    Computes the stable commutator length of the element c in the group G.

    The problem is defined by:
    - F_i = <a_i, b_i>, the free group with basis {a_i, b_i}.
    - c_i = [a_i, b_i], the commutator in F_i.
    - G is the free product of F_i for i = 1, ..., 19.
    - c is the product product_{i=1 to 19} c_i^30 in G.
    """

    # Number of free groups in the free product
    num_groups = 19

    # The exponent of each commutator
    exponent = 30

    # The stable commutator length (scl) of a single commutator [a, b] in a
    # free group F_2 = <a, b> is a known value.
    scl_of_one_commutator = 0.5

    # Step 1: Compute the scl of one term c_i^30 in its corresponding group F_i.
    # Using the homogeneity property scl(g^n) = n * scl(g).
    scl_of_one_powered_term = exponent * scl_of_one_commutator

    # Step 2: Compute the total scl of c in G.
    # The scl is additive over free products.
    total_scl = num_groups * scl_of_one_powered_term

    # Print the final equation with all the numbers.
    print("The final calculation for the stable commutator length is as follows:")
    print(f"Total scl = (Number of groups) * (Exponent) * (scl of a single commutator)")
    print(f"Total scl = {num_groups} * {exponent} * {scl_of_one_commutator} = {total_scl}")

if __name__ == "__main__":
    solve_scl_problem()