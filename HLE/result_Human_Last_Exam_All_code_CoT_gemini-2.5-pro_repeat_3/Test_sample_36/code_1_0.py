import math
from functools import reduce
from itertools import product, chain, combinations

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b) if math.gcd(a, b) != 0 else 0

def lcm_list(numbers):
    """Computes the least common multiple of a list of numbers."""
    if not numbers:
        return 1
    return reduce(lcm, numbers, 1)

def get_element_orders(group_structure):
    """
    Calculates the set of all possible element orders for a group
    given by a direct product of cyclic groups.
    group_structure: A list of orders of the cyclic components, e.g., [2, 3, 3].
    """
    divisor_sets = [get_divisors(k) for k in group_structure]
    
    element_orders = set()
    # The order of an element (g1, g2, ...) is lcm(ord(g1), ord(g2), ...)
    # where ord(gi) must be a divisor of the order of the cyclic group Ci.
    for combo in product(*divisor_sets):
        if not combo:
            continue
        order = lcm_list(combo)
        element_orders.add(order)
        
    return element_orders

def get_all_non_empty_subsets(s):
    """Returns all non-empty subsets of a set s."""
    ss = list(s)
    return chain.from_iterable(combinations(ss, r) for r in range(1, len(ss) + 1))

def main():
    """
    Calculates the number of unique sets S(rho) for Abelian groups of order 18.
    """
    # The two Abelian groups of order 18 are C_18 and C_2 x C_3 x C_3.
    # We represent them by the orders of their cyclic components.
    group_structures = [
        [18],                     # Represents C_18
        [2, 3, 3]                 # Represents C_2 x C_3 x C_3
    ]
    
    all_possible_set_orders = set()
    
    print("Analyzing Abelian groups of order 18:")
    print("-" * 40)
    
    for i, struct in enumerate(group_structures, 1):
        # The orders of characters are the same as the orders of elements.
        element_orders = get_element_orders(struct)
        
        # The order k of a set S(rho) = U_k is the lcm of the orders of the
        # characters in the representation rho. We find all possible such k.
        possible_k = set()
        for subset in get_all_non_empty_subsets(element_orders):
            k = lcm_list(subset)
            possible_k.add(k)
        
        group_name = " x ".join([f"C_{c}" for c in struct])
        print(f"Group {i}: {group_name}")
        print(f"  - Element/Character orders: {sorted(list(element_orders))}")
        print(f"  - Possible orders k for S(rho)=U_k: {sorted(list(possible_k))}")
        print("-" * 40)
        
        all_possible_set_orders.update(possible_k)
        
    print(f"Union of all possible orders k: {sorted(list(all_possible_set_orders))}")
    final_count = len(all_possible_set_orders)
    print(f"Total number of unique sets S(rho) is the count of these unique orders.")
    print(f"Final Answer: {final_count}")

if __name__ == '__main__':
    main()