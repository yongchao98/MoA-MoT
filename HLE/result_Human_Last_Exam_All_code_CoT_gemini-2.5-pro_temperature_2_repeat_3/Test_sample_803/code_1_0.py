import libgap
import sys

def find_filled_groups(q, m):
    """
    Finds and prints all nonabelian filled groups of order 2*q**m.

    A nonabelian group G of order 2*q**m (q odd prime) is filled if and only if
    the order of its abelianization G/[G,G] is 2.

    Args:
        q (int): An odd prime number.
        m (int): A natural number.
    """
    if not (q > 2 and all(q % i != 0 for i in range(2, int(q**0.5) + 1))):
        print(f"Error: q={q} is not an odd prime.", file=sys.stderr)
        return
    if not (m > 0 and isinstance(m, int)):
        print(f"Error: m={m} is not a natural number.", file=sys.stderr)
        return
        
    n = 2 * (q ** m)
    
    print(f"Searching for nonabelian filled groups of order n = 2 * {q}^{m} = {n}")
    
    try:
        num_groups = libgap.NrSmallGroups(n)
        if num_groups == 0:
             print(f"The SmallGroups library in GAP does not contain groups of order {n}.")
             return
    except Exception as e:
        print(f"Could not retrieve groups of order {n}. It might be too large for the SmallGroups library.", file=sys.stderr)
        print(f"GAP error: {e}", file=sys.stderr)
        return
        
    print(f"Found {num_groups} groups of order {n} in the library.")
    print("---")
    
    filled_groups_found = []
    
    for i in range(1, num_groups + 1):
        # Access the group from the SmallGroups library
        G = libgap.SmallGroup(n, i)
        
        # We are only interested in nonabelian groups
        if G.IsAbelian():
            continue
            
        # The condition for being a filled group of this order is |G/[G,G]| = 2.
        # G.CommutatorSubgroup() could also be used.
        derived_subgroup = G.DerivedSubgroup()
        abelianization_order = G.Order() / derived_subgroup.Order()
        
        if abelianization_order == 2:
            group_id = f"SmallGroup({n}, {i})"
            structure_desc = G.StructureDescription()
            # The final equation is the group's structure description.
            # Example: G = C_9 : C_2
            final_equation = f"G = {structure_desc}"
            
            output = f"Found filled group: {group_id}\n"
            output += f"Structure: {final_equation}\n"
            # Output each number from the structure description equation
            numbers_in_equation = [c for c in structure_desc if c.isdigit()]
            output += f"Numbers in structure: {', '.join(numbers_in_equation)}"
            
            filled_groups_found.append(output)

    if not filled_groups_found:
        print("No nonabelian filled groups found for the given q and m.")
    else:
        print("\n\n".join(filled_groups_found))
        
    print("---")


if __name__ == '__main__':
    # Example: Find nonabelian filled groups of order 18 (q=3, m=2)
    find_filled_groups(q=3, m=2)
    
    print("\n=============================================\n")
    
    # Example: Find nonabelian filled groups of order 50 (q=5, m=2)
    find_filled_groups(q=5, m=2)