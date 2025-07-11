import math
from functools import reduce
from itertools import product

def get_divisors(n):
    """Returns a list of all positive divisors of an integer n."""
    divs = set()
    # 1 and n are always divisors
    divs.add(1)
    divs.add(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b)

def get_element_orders(group_structure):
    """
    Computes the set of all possible element orders for a group G = Z_n1 x Z_n2 x ...
    'group_structure' is a list of integers [n1, n2, ...].
    """
    divisor_sets = [get_divisors(n) for n in group_structure]
    order_combinations = product(*divisor_sets)
    
    all_orders = set()
    for orders_tuple in order_combinations:
        current_lcm = reduce(lcm, orders_tuple, 1)
        all_orders.add(current_lcm)
    return all_orders

def solve_and_print():
    """
    Solves the problem for Abelian groups of order 18 and prints the steps.
    """
    # The two non-isomorphic Abelian groups of order 18
    # G1 is isomorphic to Z_18
    group1_struct = [18]
    # G2 is Z_2 x Z_3 x Z_3
    group2_struct = [2, 3, 3]

    print("The two non-isomorphic Abelian groups of order 18 are G1 = Z_18 and G2 = Z_2 x Z_3 x Z_3.")

    # Calculate orders for G1 = Z_18
    orders_g1 = get_element_orders(group1_struct)
    sorted_orders_g1 = sorted(list(orders_g1))
    print(f"\nThe set of element orders for G1 (Z_18) are the divisors of 18:")
    print(sorted_orders_g1)

    # Calculate orders for G2 = Z_2 x Z_3 x Z_3
    orders_g2 = get_element_orders(group2_struct)
    sorted_orders_g2 = sorted(list(orders_g2))
    print(f"\nThe set of element orders for G2 (Z_2 x Z_3 x Z_3) are the lcm of divisors of 2, 3, and 3:")
    print(sorted_orders_g2)
    
    # The total unique sets correspond to the union of these sets of orders.
    all_possible_orders = orders_g1.union(orders_g2)
    sorted_all_orders = sorted(list(all_possible_orders))

    print("\nThe combined set of unique orders is the union of the two sets above.")
    print(f"Union = {sorted_orders_g1} U {sorted_orders_g2}")
    print(f"Result = {sorted_all_orders}")

    final_count = len(all_possible_orders)
    print(f"\nEach of these {final_count} orders corresponds to a unique set of eigenvalues (the k-th roots of unity).")
    print("Final equation:")
    # Printing each number in the final equation
    numbers_str = ", ".join(map(str, sorted_all_orders))
    print(f"count = size of {{{numbers_str}}} = {final_count}")

solve_and_print()