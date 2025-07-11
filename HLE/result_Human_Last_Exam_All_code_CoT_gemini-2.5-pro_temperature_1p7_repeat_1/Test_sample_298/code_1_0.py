def get_group_str(free_rank, torsion_parts):
    """Formats a group into the specified string format."""
    parts = []
    if free_rank > 0:
        if free_rank == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{free_rank}")

    for t in sorted(torsion_parts):
        parts.append(f"Z/{t}Z")
    
    if not parts:
        return "0"
    
    return "+".join(parts)

def solve_cohomology():
    """
    Calculates and prints the list of cohomology groups for M(7).
    """
    # Step 1 & 2: Set up the homology groups for Conf_7(S^1)
    # H_i = (free_rank, [torsion_parts])
    homology_groups = {
        0: (1, []),
        1: (1, [6]),
        2: (0, [60, 2]),
        3: (0, [12, 6]),
        4: (0, [12, 6]),
        5: (0, [30, 2, 2]),
        6: (0, []),
    }
    
    max_dim = max(homology_groups.keys())
    
    cohomology_groups = []
    
    # Step 3: Apply Universal Coefficient Theorem
    # H^i = Free(H_i) + Torsion(H_{i-1})
    
    # H^0 = Free(H_0)
    H0_free_rank, H0_torsion = homology_groups.get(0, (0, []))
    cohomology_groups.append(get_group_str(H0_free_rank, []))

    for i in range(1, max_dim + 2):
        Hi_free_rank, Hi_torsion = homology_groups.get(i, (0, []))
        H_prev_free_rank, H_prev_torsion = homology_groups.get(i - 1, (0, []))
        
        current_cohomology_str = get_group_str(Hi_free_rank, H_prev_torsion)
        cohomology_groups.append(current_cohomology_str)
        
    # Remove trailing zeros
    last_nonzero = -1
    for i in range(len(cohomology_groups) - 1, -1, -1):
        if cohomology_groups[i] != "0":
            last_nonzero = i
            break
            
    final_list = cohomology_groups[:last_nonzero + 1]
    
    print(f"[{', '.join(final_list)}]")


solve_cohomology()