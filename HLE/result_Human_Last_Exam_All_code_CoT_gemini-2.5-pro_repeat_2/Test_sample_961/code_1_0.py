def solve_scl():
    """
    Calculates the stable commutator length based on the problem description.
    """
    # Number of free groups in the free product, G = F_1 * ... * F_19.
    num_groups = 19

    # The exponent of each commutator c_i in the product c.
    exponent = 30

    # The stable commutator length (scl) of a basic commutator [a, b]
    # in a free group F = <a, b> is a standard result, scl_F([a, b]) = 1/2.
    scl_basic_commutator = 0.5

    # The calculation follows from properties of scl:
    # scl_G(c) = sum_{i=1 to 19} scl_{F_i}(c_i^30)  (additivity in free products)
    #          = sum_{i=1 to 19} 30 * scl_{F_i}(c_i) (homogeneity)
    #          = 19 * 30 * scl_F([a,b])             (since all terms are identical)
    #          = 19 * 30 * 0.5
    final_scl = num_groups * exponent * scl_basic_commutator

    # Ensure the final result is an integer if it's a whole number for cleaner output.
    if final_scl.is_integer():
        final_scl = int(final_scl)

    print("The stable commutator length of c is calculated by the equation:")
    print(f"{num_groups} * {exponent} * {scl_basic_commutator} = {final_scl}")

solve_scl()