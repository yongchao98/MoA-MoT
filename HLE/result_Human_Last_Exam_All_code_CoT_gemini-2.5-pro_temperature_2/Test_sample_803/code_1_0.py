import math

def describe_nonabelian_filled_groups(q, m):
    """
    Identifies and describes the nonabelian filled groups of order 2 * q^m.

    This function provides a detailed classification based on established theorems
    in finite group theory. According to these results, the nonabelian filled
    groups of this order are the "generalized dihedral groups," constructed
    from the abelian groups of order q^m.

    Args:
        q (int): An odd prime number.
        m (int): A natural number (m >= 1).
    """

    # Helper function to check for primality
    def is_prime(n):
        if n <= 1: return False
        if n <= 3: return True
        if n % 2 == 0 or n % 3 == 0: return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    # Helper function to generate integer partitions
    def find_partitions(n):
        partitions_list = []
        def generate(target, max_part, current_part):
            if target == 0:
                partitions_list.append(current_part)
                return
            for i in range(min(target, max_part), 0, -1):
                generate(target - i, i, current_part + [i])
        generate(n, n, [])
        return partitions_list

    # Input validation
    if not isinstance(q, int) or not isinstance(m, int):
        print("Error: q and m must be integers.")
        return
    if q < 3 or q % 2 == 0 or not is_prime(q):
        print(f"Error: The provided value q={q} is not an odd prime number.")
        return
    if m < 1:
        print(f"Error: The provided value m={m} must be a natural number (m >= 1).")
        return
        
    group_order_val = 2 * (q**m)
    h_order_val = q**m
    
    partitions = find_partitions(m)
    num_groups = len(partitions)
    
    print(f"For q={q} and m={m}:")
    print(f"The nonabelian filled groups of order 2 * {q}^{m} = {group_order_val} are generalized dihedral groups.")
    print(f"They are constructed from the abelian groups H of order {q}^{m} = {h_order_val}.")
    print("These groups G are semidirect products of the form G = H ⋊ Z₂, where the non-trivial")
    print("element of Z₂ acts by sending every element h ∈ H to its inverse h⁻¹.\n")
    print(f"There are {num_groups} such non-isomorphic group(s), corresponding to the {num_groups} partition(s) of m = {m}.\n")
    print("-" * 60)

    for i, p in enumerate(partitions):
        print(f"Group {i+1} (based on partition {p}):")
        
        # Build string representations for the abelian group H
        h_symbolic = " x ".join([f"Z_{{q^{part}}}" for part in p])
        h_numeric = " x ".join([f"Z_{{{q**part}}}" for part in p])
        
        # Identify the standard dihedral group case (when H is cyclic)
        if len(p) == 1:
            h_string = f"H = Z_{{q^{m}}} = Z_{{{h_order_val}}}"
            g_string = f"G = Dih(H) = D_{{{2 * h_order_val}}} (the standard dihedral group of order {2*h_order_val})"
        else:
            h_string = f"H = {h_symbolic} = {h_numeric}"
            g_string = f"G = Dih(H) = Dih({h_numeric})"
        
        print(f"  Abelian Group Component: {h_string}")
        print(f"  Nonabelian Filled Group: {g_string}\n")


# --- Main execution block ---
# The user can modify the values of q and m here to explore different cases.

# Example 1: The smallest case for a non-dihedral group of this type
print("Demonstration for q=3, m=2 (order 18):")
describe_nonabelian_filled_groups(q=3, m=2)

print("\n" + "="*60 + "\n")

# Example 2: A case with more partitions
print("Demonstration for q=5, m=3 (order 250):")
describe_nonabelian_filled_groups(q=5, m=3)
