import math

def is_p_group_filled(p_group_name):
    """
    Checks if a given p-group is filled based on established classifications.

    A p-group is filled if and only if:
    1. It has exponent p (every non-identity element has order p).
    2. It is the quaternion group Q8.

    Args:
        p_group_name (str): A string representing the p-group.
                            Examples: 'C5', 'C2 x C2 x C2', 'Q8', 'D8', 'C9'.

    Returns:
        bool: True if the p-group is filled, False otherwise.
    """
    # Case 1: Check for the special case of Q8
    if p_group_name == 'Q8':
        return True

    # A small database of properties for common p-groups
    # (name, p, exponent)
    p_group_db = {
        # Elementary abelian groups (exponent p) -> Filled
        'C2': (2, 2),
        'C3': (3, 3),
        'C5': (5, 5),
        'C7': (7, 7),
        'C2 x C2': (2, 2),
        'C3 x C3': (3, 3),
        # Quaternion group -> Filled
        'Q8': (2, 4), # Special case handled above, but included for completeness
        # Groups with exponent > p -> Not Filled
        'C4': (2, 4),
        'C8': (2, 8),
        'C9': (3, 9),
        'C25': (5, 25),
        'C2 x C4': (2, 4),
        'D8': (2, 4), # Dihedral group of order 8
    }

    if p_group_name not in p_group_db:
        print(f"Warning: Information for p-group '{p_group_name}' is not available in the database.")
        print("Assuming it's not a special case and checking its exponent if possible.")
        # A simple parser for C_p^k style names
        if p_group_name.startswith('C'):
            try:
                order = int(p_group_name[1:])
                is_power_of_prime = False
                if order > 1:
                    # Check if order is p^k
                    p = None
                    temp_order = order
                    for i in range(2, int(math.sqrt(order)) + 2):
                        if temp_order % i == 0:
                            p = i
                            while temp_order % p == 0:
                                temp_order //= p
                            if temp_order == 1:
                                is_power_of_prime = True
                            break
                    if temp_order > 1 and p is None: # order is prime
                         p = temp_order
                         is_power_of_prime = True

                if is_power_of_prime:
                    # For C_{p^k}, exponent is p^k
                    # It's filled iff exponent is p, i.e., k=1
                    return order == p
            except ValueError:
                pass # Fallback to False
        return False

    p, exponent = p_group_db[p_group_name]

    # Case 2: Check if the exponent is p
    return exponent == p

def check_finite_filled_nilpotent_group(sylow_subgroups):
    """
    Determines if a finite nilpotent group is a filled group.

    A finite nilpotent group is the direct product of its Sylow p-subgroups.
    It is filled if and only if all of its Sylow p-subgroups are filled.

    Args:
        sylow_subgroups (list): A list of strings, where each string is the
                                name of a Sylow p-subgroup.
                                e.g., ['C3', 'Q8'] for the group C3 x Q8.
    """
    print(f"Analyzing the nilpotent group: {' x '.join(sylow_subgroups)}")
    
    all_sylows_filled = True
    for p_group in sylow_subgroups:
        is_filled = is_p_group_filled(p_group)
        print(f"- Checking Sylow subgroup '{p_group}': {'Filled' if is_filled else 'Not Filled'}")
        if not is_filled:
            all_sylows_filled = False

    print("\n--- Conclusion ---")
    if all_sylows_filled:
        print("The group is a FILLED group.")
    else:
        print("The group is NOT a filled group.")
    print("A finite nilpotent group is filled if and only if each of its Sylow p-subgroups P")
    print("has exponent p (meaning every element x in P satisfies x^p = e) OR is the quaternion group Q8.")

if __name__ == '__main__':
    # Example 1: A group whose Sylow subgroups are C5 and Q8.
    # C5 is a 5-group with exponent 5 (filled).
    # Q8 is a 2-group that is a known special case (filled).
    # Since all Sylow subgroups are filled, the direct product is filled.
    group1 = ['C5', 'Q8']
    check_finite_filled_nilpotent_group(group1)

    print("\n" + "="*40 + "\n")

    # Example 2: A group whose Sylow subgroups are C3 and D8.
    # C3 is a 3-group with exponent 3 (filled).
    # D8 is a 2-group with exponent 4, which is not 2, and it's not Q8 (not filled).
    # Since one Sylow subgroup is not filled, the group is not filled.
    group2 = ['C3', 'D8']
    check_finite_filled_nilpotent_group(group2)
    
    print("\n" + "="*40 + "\n")

    # Example 3: The abelian group C12 = C3 x C4
    # C3 is filled.
    # C4 is a 2-group with exponent 4, which is not 2 (not filled).
    # The group is not filled.
    group3 = ['C3', 'C4']
    check_finite_filled_nilpotent_group(group3)