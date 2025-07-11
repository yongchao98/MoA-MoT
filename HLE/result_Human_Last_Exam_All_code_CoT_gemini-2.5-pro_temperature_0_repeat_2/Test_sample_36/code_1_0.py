import math
from functools import reduce
from itertools import product, combinations

def get_lcm(a, b):
    """Computes the least common multiple of two integers."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def get_divisors(n):
    """Computes all divisors of an integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return divs

def get_element_orders(group_factors):
    """
    Calculates the set of all possible element orders for a group
    given as a direct product of cyclic groups (e.g., Z_6 x Z_3 -> [6, 3]).
    """
    if not group_factors:
        return {1}
    
    # Get orders of elements in each cyclic component (which are the divisors)
    list_of_divisor_sets = [get_divisors(n) for n in group_factors]
    
    element_orders = set()
    
    # The order of an element (g1, g2, ...) is lcm(ord(g1), ord(g2), ...)
    for ord_combination in product(*list_of_divisor_sets):
        order = reduce(get_lcm, ord_combination, 1)
        element_orders.add(order)
        
    return element_orders

def main():
    """
    Main function to solve the problem.
    """
    # Step 1: Define the non-isomorphic Abelian groups of order 18
    # by their cyclic factors. 18 = 2 * 9.
    # Group 1: Z_18 (isomorphic to Z_2 x Z_9)
    # Group 2: Z_6 x Z_3 (isomorphic to Z_2 x Z_3 x Z_3)
    groups = {
        "Z_18": [18],
        "Z_6 x Z_3": [6, 3]
    }
    
    print("Finding the number of unique sets S(rho) for Abelian groups of order 18.\n")

    all_possible_set_indices = set()
    
    # Process each group
    for name, factors in groups.items():
        print(f"--- Analyzing Group: {name} ---")
        
        # Step 2: Find all possible orders of elements in the group.
        element_orders = get_element_orders(factors)
        print(f"The set of element orders is: {sorted(list(element_orders))}")
        
        # Step 3: Find all possible lcm's of non-empty subsets of these orders.
        # This gives the indices 'd' for the possible sets U_d.
        possible_lcm_values = set()
        order_list = list(element_orders)
        for r in range(1, len(order_list) + 1):
            for subset in combinations(order_list, r):
                current_lcm = reduce(get_lcm, subset, 1)
                possible_lcm_values.add(current_lcm)
        
        print(f"The possible indices 'd' for the sets S(rho) = U_d are: {sorted(list(possible_lcm_values))}\n")
        
        # Add these to the total collection of possible sets.
        all_possible_set_indices.update(possible_lcm_values)

    # Step 4: Combine results and find the total count.
    print("--- Final Result ---")
    print("The complete set of unique indices 'd' from all groups of order 18 is the union of the sets found above.")
    final_indices = sorted(list(all_possible_set_indices))
    print(f"Union of indices = {final_indices}")
    
    final_count = len(all_possible_set_indices)
    print(f"\nThe total number of unique sets S(rho) is the size of this combined set of indices.")
    print(f"Total unique sets = {final_count}")
    
    print(f"<<<{final_count}>>>")

if __name__ == "__main__":
    main()