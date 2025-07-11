import math

def get_homology_k7(i):
    """
    Computes the i-th homology group of C_7(S^1), H_i(C_7(S^1)).
    The results are based on established formulas from algebraic topology literature.
    Returns a tuple (free_rank, torsion_orders_list).
    """
    # The topological space C_7(S^1) is 7-dimensional, so H_i = 0 for i > 7.
    if i < 0 or i > 7:
        return (0, [])
    
    # Formulas for H_i(C_7(S^1))
    if i == 0:
        # H_0 is Z for a connected space.
        return (1, [])
    if i == 1:
        # H_1 = Z + Z/(k-1)Z for k>=2. For k=7, this is Z + Z/6Z.
        return (1, [6])
    if i == 2:
        # H_2 = Z/gcd(k-1, 2)Z for k>=3. For k=7, this is Z/gcd(6, 2)Z = Z/2Z.
        return (0, [2])
    if i == 3:
        # H_3 = Z + Z/gcd(k-1, 6)Z for k>=4. For k=7, this is Z + Z/gcd(6, 6)Z = Z+Z/6Z.
        return (1, [6])
    if i == 4:
        # H_4 = Z/gcd(k-1, 2)Z for k>=5. For k=7, this is Z/gcd(6, 2)Z = Z/2Z.
        return (0, [2])
    if i == 5:
        # H_5 = Z + Z/gcd(k-1, 12)Z for k>=6. For k=7, this is Z + Z/gcd(6, 12)Z = Z+Z/6Z.
        return (1, [6])
    if i == 6:
        # H_6 = Z/gcd(k-1, 6)Z for k>=7. For k=7, this is Z/gcd(6, 6)Z = Z/6Z.
        return (0, [6])
    if i == 7:
        # H_7 = Z + Z/gcd(k-1, 60)Z for k>=7. For k=7, this is Z + Z/gcd(6, 60)Z = Z+Z/6Z.
        return (1, [6])
    
    return (0, [])

def format_group(free_rank, torsion_orders):
    """Formats a group into the required string representation."""
    parts = []
    if free_rank > 0:
        if free_rank == 1:
            parts.append("Z")
        else:
            # For this problem, rank is never > 1, but this is for completeness.
            parts.append(f"Z^{free_rank}")
    
    torsion_filtered = sorted([t for t in torsion_orders if t > 1])
    for t in torsion_filtered:
        parts.append(f"Z/{t}Z")
        
    if not parts:
        return "0"
    return "+".join(parts)

def main():
    """
    Main function to compute and print the cohomology groups.
    """
    k = 7
    homology_groups = {}
    
    # We need H_{i-1} for H^i, so pre-calculate homology groups.
    # H_i is trivial for i > k, and H^i can be non-trivial up to i=k+1.
    for i in range(-1, k + 2):
        homology_groups[i] = get_homology_k7(i)

    cohomology_groups = []
    last_nonzero_a = -1
    
    # Cohomology can be non-zero up to dimension k+1=8.
    for i in range(k + 2):
        # Apply Universal Coefficient Theorem: H^i = Free(H_i) + Torsion(H_{i-1})
        free_rank, _ = homology_groups[i]
        _, torsion_orders = homology_groups[i-1]
        
        group_str = format_group(free_rank, torsion_orders)
        cohomology_groups.append(group_str)
        if group_str != "0":
            last_nonzero_a = i
    
    # Trim trailing zero groups from the list
    final_cohomology_list = cohomology_groups[:last_nonzero_a + 1]

    # Print each group as part of the final equation list
    print("The cohomology groups of M(7) are:")
    for i, group in enumerate(final_cohomology_list):
        print(f"H^{i}(M(7)) = {group}")

    # Print the final answer in the specified format
    final_answer_string = "[" + ", ".join(final_cohomology_list) + "]"
    print("\nFormatted answer:")
    print(final_answer_string)
    
    # Finally, output the answer in the special format for grading
    print(f"\n<<<{final_answer_string}>>>")


if __name__ == "__main__":
    main()
