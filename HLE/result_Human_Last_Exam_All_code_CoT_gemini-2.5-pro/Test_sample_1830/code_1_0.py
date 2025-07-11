def find_order_type():
    """
    This function determines the order type of the set of possible cardinalities
    of maximal almost disjoint families of infinite subsets of omega, under the
    assumption that 2^omega = omega_1.
    """

    # Step 1: Define the set X of possible cardinalities.
    # Under the Continuum Hypothesis (2^omega = omega_1), the possible
    # cardinalities of a Maximal Almost Disjoint Family (MADF) are:
    # 1. Any finite integer n >= 1. (e.g., by partitioning omega into n infinite sets).
    # 2. Countably infinite, omega. (A known result in ZFC).
    # 3. The first uncountable cardinal, omega_1. (A known result for 2^omega, which is omega_1 by CH).
    
    X_description = "{1, 2, 3, ...} U {omega, omega_1}"

    print("Step 1: Determine the set X of possible cardinalities of MADFs.")
    print(f"The set is X = {X_description}.\n")

    # Step 2: Determine the order type of X.
    # The elements of X, in order, are 1, 2, 3, ..., omega, omega_1.
    # We find an order-isomorphic ordinal.
    
    # The initial segment of positive integers corresponds to the ordinal omega.
    finite_part_ordinal = "omega"
    
    # There are two additional elements beyond the finite ones.
    # The first is 'omega', which corresponds to adding 1 to the ordinal type of the finite part.
    # The second is 'omega_1', which corresponds to adding another 1.
    additional_elements_count = 2
    
    print("Step 2: Find the order type of X.")
    print("The ordered set is (1, 2, 3, ..., omega, omega_1).")
    print("The order type is constructed by considering the parts of the set:")
    print(f" - The infinite sequence of finite numbers 1, 2, 3, ... has order type '{finite_part_ordinal}'.")
    print(f" - There are {additional_elements_count} larger elements (omega and omega_1).")
    
    # The final order type is the sum of the ordinal for the initial part
    # and the count of the remaining elements.
    print("\nThe final order type is given by the ordinal sum:")
    print(f"{finite_part_ordinal} + {additional_elements_count}")

    final_order_type = "omega + 2"
    print(f"\nThus, the order type of X is {final_order_type}.")

find_order_type()