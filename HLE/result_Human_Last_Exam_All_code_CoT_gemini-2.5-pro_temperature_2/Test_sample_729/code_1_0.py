import collections

def get_q128_elements():
    """Generates all 128 elements of the group Q_128."""
    elements = set()
    # Elements of the form x^i
    for i in range(64):
        elements.add(('x', i))
    # Elements of the form x^i y
    for i in range(64):
        elements.add(('xy', i))
    return elements

def element_power(g, k):
    """
    Computes g^k for an element g in Q_128.
    The group is Q_128 = <x, y | x^64=1, x^32=y^2, y^-1xy=x^-1>.
    n is such that |G| = 4n, so n=32. The cyclic part <x> has order 2n=64.
    
    Args:
        g (tuple): An element represented as ('x', i) or ('xy', i).
        k (int): The exponent.
        
    Returns:
        tuple: The element g^k in the same representation.
    """
    g_type, i = g
    cyclic_order = 64
    n_val = 32

    if g_type == 'x':
        # (x^i)^k = x^(i*k)
        return ('x', (i * k) % cyclic_order)
    elif g_type == 'xy':
        # The element x^i y has order 4. (x^i y)^2 = x^32.
        # The result depends on k mod 4.
        k_mod_4 = k % 4
        if k_mod_4 == 0:
            return ('x', 0)  # Identity element x^0
        elif k_mod_4 == 1:
            return ('xy', i)  # (x^i y)^1 = x^i y
        elif k_mod_4 == 2:
            return ('x', n_val)  # (x^i y)^2 = x^32
        else:  # k_mod_4 == 3
            # (x^i y)^3 = (x^i y)^2 * (x^i y) = x^32 * x^i y = x^(i+32) y
            return ('xy', (i + n_val) % cyclic_order)

def find_power_subgroups():
    """
    Finds all unique power subgroups of the generalized quaternion group Q_128.
    """
    group_order = 128
    all_elements = get_q128_elements()
    
    unique_power_subgroups = set()
    
    # We check k from 0 to group_order-1, which is sufficient to find all distinct power subgroups.
    for k in range(group_order):
        # G^k = {g^k for all g in G}
        power_set = frozenset(element_power(g, k) for g in all_elements)
        unique_power_subgroups.add(power_set)

    print("Found the following unique power subgroups:")
    # Sort by size for clearer output
    sorted_subgroups = sorted(list(unique_power_subgroups), key=len, reverse=True)
    
    for i, subgroup in enumerate(sorted_subgroups):
        desc = ""
        if len(subgroup) == 128:
            desc = " (The group Q_128 itself)"
        elif len(subgroup) == 1:
            desc = " (The trivial subgroup {1})"
        else:
            is_cyclic = all(g[0] == 'x' for g in subgroup)
            if is_cyclic:
                desc = f" (Cyclic subgroup of order {len(subgroup)})"

        print(f"- Subgroup {i+1}: Found a power subgroup of size {len(subgroup)}{desc}.")

    print(f"\nIn total, there are {len(unique_power_subgroups)} distinct power subgroups in the generalized quaternion group of size 128.")

if __name__ == '__main__':
    find_power_subgroups()
