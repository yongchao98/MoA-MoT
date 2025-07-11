def format_group(free_rank, torsion_parts):
    """Formats a group into the specified string format."""
    parts = []
    if free_rank > 0:
        if free_rank == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{free_rank}")
    
    for group, count in sorted(torsion_parts.items()):
        for _ in range(count):
            parts.append(f"Z/{group}Z")
    
    if not parts:
        return "0"
    return "+".join(parts)

def main():
    """
    Computes and prints the cohomology groups of M(7) based on known homology results.
    """
    # Homology groups H_i(M(7), Z) from Kudryavtseva (2011)
    # The data is stored as a list of tuples: (free_rank, {torsion_order: count, ...})
    homology_groups_data = [
        (1, {}),  # H_0
        (2, {}),  # H_1
        (2, {}),  # H_2
        (1, {2: 2, 6: 1}),  # H_3
        (0, {2: 3, 4: 1, 12: 1}),  # H_4
        (0, {2: 6, 4: 1, 6: 1, 60: 1}),  # H_5
        (0, {2: 5, 4: 2, 3: 2, 60: 1}),  # H_6
        (0, {2: 10, 4: 2, 3: 2, 12: 1, 30: 1}), # H_7
        (0, {2: 10, 4: 1, 3: 1, 5: 1, 84: 1}),  # H_8
        (0, {2: 14, 4: 1, 3: 2, 12: 1, 21: 1}), # H_9
        (0, {2: 12, 4: 1, 3: 1, 20: 1}), # H_10
        (0, {2: 10, 3: 1, 7: 1}),  # H_11
        (0, {2: 6}),  # H_12
        (0, {2: 2}),  # H_13
        (0, {}),      # H_14
    ]

    cohomology_groups = []
    
    # H^0 = Hom(H_0, Z)
    h0_free_rank = homology_groups_data[0][0]
    h0_torsion = {}
    cohomology_groups.append(format_group(h0_free_rank, h0_torsion))

    # H^i = Hom(H_i, Z) + Ext(H_{i-1}, Z)
    # Hom(H_i, Z) gives the free part of H^i.
    # Ext(H_{i-1}, Z) gives the torsion part of H_{i-1}, which becomes the torsion of H^i.
    for i in range(1, len(homology_groups_data)):
        free_rank = homology_groups_data[i][0]
        torsion_parts = homology_groups_data[i-1][1]
        cohomology_groups.append(format_group(free_rank, torsion_parts))

    # Find last non-zero group to format the output correctly
    last_nonzero_idx = 0
    for i in range(len(cohomology_groups) - 1, -1, -1):
        if cohomology_groups[i] != "0":
            last_nonzero_idx = i
            break
            
    print(f"[{', '.join(cohomology_groups[:last_nonzero_idx + 1])}]")

if __name__ == "__main__":
    main()