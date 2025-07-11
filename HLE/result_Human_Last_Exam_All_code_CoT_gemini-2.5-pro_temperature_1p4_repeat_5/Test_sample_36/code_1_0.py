import math
from functools import reduce
from itertools import product

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return abs(a)

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def list_lcm(numbers):
    """Computes the least common multiple of a list of numbers."""
    if not numbers:
        return 1  # By convention, lcm of an empty set is 1.
    return reduce(lcm, numbers, 1)

def get_divisors(n):
    """Returns a list of all divisors of a positive integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def get_element_orders_for_group(cyclic_factors):
    """
    Calculates the set of all possible element orders for an Abelian group
    given by its representation as a product of cyclic groups, e.g., [2, 3, 3] for Z_2 x Z_3 x Z_3.
    """
    # The order of an element (g1, g2, ...) is lcm(ord(g1), ord(g2), ...).
    # The possible orders for an element in Z_n are the divisors of n.
    list_of_divisors = [get_divisors(f) for f in cyclic_factors]
    
    element_orders = set()
    # The cartesian product gives all combinations of divisor choices.
    for combo in product(*list_of_divisors):
        order = list_lcm(combo)
        element_orders.add(order)
        
    return element_orders

def solve_for_order_18():
    """
    Main function to solve the problem for Abelian groups of order 18.
    """
    order = 18
    # From the fundamental theorem of finite Abelian groups, the groups of order 18 ($2 \times 3^2$) are:
    # 1. G1 = Z_18 (isomorphic to Z_2 x Z_9)
    # 2. G2 = Z_2 x Z_3 x Z_3
    groups = {
        "Z_18": [18],
        "Z_2 x Z_3 x Z_3": [2, 3, 3]
    }
    
    print(f"This script calculates the number of unique eigenvalue sets S(rho) for representations of Abelian groups of order {order}.\n")
    print("Step 1: The two Abelian groups of order 18 are Z_18 and Z_2 x Z_3 x Z_3.\n")
    print("Step 2: For each group G, we find the set of possible orders k for the unique eigenvalue sets U_k.")
    print("This set is the set of all element orders in G, denoted E(G).\n")

    all_k_values = set()
    
    for name, factors in groups.items():
        print(f"--- Analyzing Group G = {name} ---")
        
        # Calculate the set of element orders for the group.
        possible_k_for_group = get_element_orders_for_group(factors)
        
        print(f"The set of element orders E(G) is: {sorted(list(possible_k_for_group))}")
        print("These are the possible orders k for the sets of eigenvalues U_k for this group.")
        all_k_values.update(possible_k_for_group)
        print("-" * 35 + "\n")

    print("--- Step 3: Combine Results ---")
    final_set_of_orders = sorted(list(all_k_values))
    final_count = len(final_set_of_orders)
    
    print("The total set of possible orders k is the union of the sets from all groups.")
    print(f"Union of orders k = {final_set_of_orders}")
    
    print(f"\nEach unique order k corresponds to a unique set of eigenvalues U_k.")
    print(f"The total number of unique sets is the size of this combined set of orders.")
    print(f"Final count = {final_count}")

# Run the solver
solve_for_order_18()