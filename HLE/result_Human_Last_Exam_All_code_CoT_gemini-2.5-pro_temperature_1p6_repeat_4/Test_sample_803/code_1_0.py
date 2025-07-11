# The 'sympy' library is required to run this code.
# You can install it using your terminal: pip install sympy

import sympy
from sympy.combinatorics.named_groups import DihedralGroup

def find_and_verify_filled_group(q, m):
    """
    This function identifies the nonabelian filled group of order 2*q**m
    and verifies its properties using the sympy library.

    A finite group G is filled if it does not contain a subgroup isomorphic
    to C_p x C_p for any odd prime p. For a group of order 2*q**m, this
    means its Sylow q-subgroup must be cyclic.

    The only nonabelian group structure satisfying this is the Dihedral group.
    """
    order = 2 * (q**m)
    group_name = f"D_{order}"
    
    # As requested, output the final result in an equation format
    print(f"The nonabelian filled group of order 2 * {q}^{m} = {order} is the Dihedral Group {group_name}.")
    print("\n--- Verification using sympy ---")

    # In sympy, DihedralGroup(n) corresponds to the dihedral group of order 2n.
    # Our group D_{2*q**m} is DihedralGroup(q**m) in sympy.
    n = q**m
    try:
        G = DihedralGroup(n)
    except Exception as e:
        print(f"Error: Could not construct DihedralGroup({n}) in sympy. {e}")
        return

    # Property 1: Check Order
    print(f"1. Checking order of {group_name}:")
    print(f"   Expected: {order}, Actual: {G.order()}")
    if G.order() != order:
        print("   Verification FAILED: Order mismatch.")
        return

    # Property 2: Check if Nonabelian
    print(f"\n2. Checking if {group_name} is nonabelian:")
    is_abelian = G.is_abelian
    print(f"   Is the group abelian? {is_abelian}")
    if is_abelian:
        print("   Verification FAILED: Group is abelian.")
        return

    # Property 3: Check the 'filled' condition (Sylow q-subgroup is cyclic)
    print(f"\n3. Checking the 'filled' condition for {group_name}:")
    print(f"   This requires the Sylow {q}-subgroup to be cyclic.")
    sylow_q_subgroup = G.sylow_subgroup(q)
    is_sylow_q_cyclic = sylow_q_subgroup.is_cyclic
    print(f"   Is the Sylow {q}-subgroup (order {sylow_q_subgroup.order()}) cyclic? {is_sylow_q_cyclic}")
    if not is_sylow_q_cyclic:
        print("   Verification FAILED: Sylow q-subgroup is not cyclic.")
        return

    print("\nVerification successful. All properties match the theoretical result.")


if __name__ == '__main__':
    # --- Main Execution ---
    # Let's solve for a specific case: an odd prime q=3 and natural number m=2.
    # The order of the group is 2 * 3^2 = 18.
    q_val = 3
    m_val = 2
    find_and_verify_filled_group(q_val, m_val)