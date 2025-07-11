import math

def solve_scl():
    """
    Computes the stable commutator length of the given element.
    """
    # Number of free groups in the free product
    num_groups = 19

    # Exponent of the commutators
    exponent = 30

    # The stable commutator length of a basic commutator c_i = [a_i, b_i]
    # in a free group F_i = <a_i, b_i> is 1/2.
    scl_single_commutator = 0.5

    print("Step 1: The stable commutator length (scl) is additive over free products.")
    print(f"scl(c) = SUM_{{i=1 to {num_groups}}} scl(c_i^{{{exponent}}})")
    print("")

    print("Step 2: The scl is homogeneous, meaning scl(g^n) = n * scl(g).")
    scl_of_powered_commutator = exponent * scl_single_commutator
    print(f"scl(c_i^{{{exponent}}}) = {exponent} * scl(c_i) = {exponent} * {scl_single_commutator} = {scl_of_powered_commutator}")
    print("")

    print("Step 3: The total scl is the sum of 19 identical terms.")
    total_scl = num_groups * scl_of_powered_commutator
    print(f"Total scl = {num_groups} * {scl_of_powered_commutator}")
    print("")

    print("Final equation:")
    # Ensure integers are used for the final print representation
    print(f"{num_groups} * {int(scl_of_powered_commutator)} = {int(total_scl)}")

if __name__ == "__main__":
    solve_scl()
