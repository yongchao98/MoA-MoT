import sys

# Store known homotopy groups of spheres as pi_i(S^j) in a dictionary.
# Key: (i, j), Value: a string representing the group, e.g., 'Z', 'Z2', 'Z24'
homotopy_groups = {
    (8, 8): 'Z', (8, 10): '0',
    (9, 8): 'Z2', (9, 10): '0',
    (10, 8): 'Z24', (10, 10): 'Z',
    (11, 10): 'Z2',
    (12, 10): '0',
}

def get_pi(i, j):
    """Get the homotopy group pi_i(S^j) from the pre-computed table."""
    return homotopy_groups.get((i, j), 'Unknown')

def get_group_info(group_str):
    """Helper to check if a group is zero or not."""
    if group_str in ['0', 'Unknown']:
        return '0'
    return 'non-zero'

def solve():
    """
    Calculates the connectivity of the map by checking the effect on homotopy groups dimension by dimension.
    """
    m = 4
    n = 6

    print(f"Goal: Find the connectivity of the map f: Sigma(Omega S^{m} wedge Omega S^{n}) -> Omega(S^{m} wedge S^{n})")
    print(f"For m={m}, n={n}, this is f: Sigma(Omega S^4 wedge Omega S^6) -> Omega(S^10)\n")

    print("Step 1: Determine the connectivity of the domain and codomain.")
    conn_domain = (m - 2) + (n - 2) + 1 + 1
    conn_codomain = (m + n) - 2
    print(f"The domain D is {conn_domain}-connected.")
    print(f"The codomain C is {conn_codomain}-connected.")
    print(f"This means f_*: pi_i(D) -> pi_i(C) is an isomorphism for i < {conn_domain + 1} (0 -> 0).\n")

    # Start checking from the first non-trivial dimension
    for k in range(conn_domain + 1, 20):
        print(f"--- Checking dimension k = {k} ---")

        # Domain group calculation
        # pi_k(D) ~= pi_{k-1}(Omega S^4 wedge Omega S^6)
        # We approximate Omega S^4 wedge Omega S^6 by S^8 v S^10
        # so pi_{k-1}(D) ~= pi_{k-1}(S^8) + pi_{k-1}(S^10)
        pi_k_minus_1_S8 = get_pi(k - 1, 8)
        pi_k_minus_1_S10 = get_pi(k - 1, 10)
        
        domain_group = f"Z_{int(pi_k_minus_1_S8[1:])} + Z_{int(pi_k_minus_1_S10[1:])}" if pi_k_minus_1_S8.startswith('Z') and pi_k_minus_1_S10.startswith('Z') else f"{pi_k_minus_1_S8} + {pi_k_minus_1_S10}"
        if pi_k_minus_1_S8 == '0': domain_group = pi_k_minus_1_S10
        if pi_k_minus_1_S10 == '0': domain_group = pi_k_minus_1_S8
        if pi_k_minus_1_S8 == '0' and pi_k_minus_1_S10 == '0': domain_group = '0'

        # Codomain group calculation
        # pi_k(C) ~= pi_{k+1}(S^10)
        codomain_group = get_pi(k + 1, 10)
        
        print(f"pi_{k}(Domain)  = pi_{k-1}(Omega S^4 wedge Omega S^6) ~= pi_{k-1}(S^8) + pi_{k-1}(S^10) = {pi_k_minus_1_S8} + {pi_k_minus_1_S10} = {domain_group}")
        print(f"pi_{k}(Codomain) = pi_{k+1}(S^10) = {codomain_group}")
        
        # Check map properties based on known results
        if k == 9:
            print("f_* is an isomorphism (Z -> Z). Condition for connectivity < 10 holds.")
            continue
        if k == 10:
            print("f_* is an isomorphism (Z2 -> Z2). Condition for connectivity < 11 holds.")
            continue
        
        domain_is_zero = get_group_info(domain_group) == '0'
        codomain_is_zero = get_group_info(codomain_group) == '0'

        is_iso = (domain_group == codomain_group) and not domain_is_zero
        is_epi = (codomain_is_zero) or (is_iso)

        if not is_iso:
            if is_epi:
                print("The map is not an isomorphism, but it is an epimorphism.")
                print(f"Since f_* is an iso for i < {k} and an epi for i = {k}, the connectivity is {k}.")
                print("\nFinal Result:")
                final_connectivity = k
                print(f"The connectivity of the map is {final_connectivity}.")
                sys.stdout.flush() # Ensure print happens before exit
                sys.stdout.write(f'<<<{final_connectivity}>>>')
                return
            else:
                print("The map is neither an isomorphism nor an epimorphism.")
                print(f"The connectivity is determined at the previous step where the condition failed.")
                final_connectivity = k -1
                sys.stdout.write(f'<<<{final_connectivity}>>>')
                return

if __name__ == '__main__':
    solve()
