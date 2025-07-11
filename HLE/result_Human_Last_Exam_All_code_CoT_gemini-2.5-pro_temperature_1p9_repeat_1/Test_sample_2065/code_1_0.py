import math

def calculate_circulons():
    """
    Calculates the number of circulon-type defects for G=SO(3)
    in d = 1 through 6 spatial dimensions.
    """

    # Cardinalities of the homotopy groups of SO(3)
    # pi_k(SO(3)) -> size
    # pi_0: group is trivial (1 element)
    # pi_1: group is Z_2 (2 elements)
    # pi_2: group is trivial (1 element)
    # pi_3: group is Z (infinite elements)
    # pi_4: group is Z_2 (2 elements)
    # pi_5: group is Z_2 (2 elements)
    homotopy_group_sizes = {
        0: 1,
        1: 2,
        2: 1,
        3: math.inf,
        4: 2,
        5: 2
    }

    # Size of pi_1(SO(3)), which is the first part of our product formula
    pi_1_size = homotopy_group_sizes[1]
    
    print("Calculating the number of circulons for gauge group G=SO(3)...")
    print("-" * 60)

    for d in range(1, 7):
        if d == 1:
            # A circle cannot be embedded in 1D space.
            result = 0
            print(f"For d={d}: A circle cannot be embedded in R^1. Number of circulons = {result}")
            continue

        # The classification group is pi_1(SO(3)) x pi_{d-1}(SO(3))
        # The number of types is the cardinality of this group.
        k = d - 1
        pi_k_size = homotopy_group_sizes[k]

        # Use strings for formatting infinity
        pi_1_size_str = str(pi_1_size)
        pi_k_size_str = 'infinity' if pi_k_size == math.inf else str(pi_k_size)
        
        equation = f"|pi_1(SO(3))| x |pi_{k}(SO(3))| = {pi_1_size_str} x {pi_k_size_str}"

        if pi_k_size == math.inf:
            result = 'infinity'
        else:
            result = pi_1_size * pi_k_size

        print(f"For d={d}: Number of circulons = {equation} = {result}")

if __name__ == '__main__':
    calculate_circulons()
