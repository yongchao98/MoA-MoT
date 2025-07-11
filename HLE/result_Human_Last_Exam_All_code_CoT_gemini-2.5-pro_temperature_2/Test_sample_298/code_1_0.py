def format_group(group_list):
    """Formats a list of integers into the required Z/nZ notation."""
    if not group_list:
        return "0"
    
    parts = []
    if 'free' in group_list:
        free_rank = group_list.count('free')
        if free_rank == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{free_rank}")
        group_list = [item for item in group_list if item != 'free']

    torsion_parts = []
    if group_list:
        for n in sorted(group_list):
            torsion_parts.append(f"Z/{n}Z")

    return "+".join(parts + torsion_parts)


def compute_cohomology_list(k):
    """
    Computes and prints the list of cohomology groups for M(k).
    For k=7, assumes H_i(B_7) = 0 for i >= 7.
    """
    if k != 7:
        print("This script is specifically configured for k=7.")
        return

    # H_i(B_7, Z) based on published results, truncated at i=6.
    # We represent Z as ['free'] and Z/nZ as [n].
    # Direct sum is represented by extending the list.
    homology = {
        0: ['free'],
        1: ['free'],
        2: [2],
        3: [12],
        4: [6],
        5: [12],
        6: [2, 2]
    }

    cohomology = {}

    # H^0
    cohomology[0] = homology[0]  # Z

    # H^1
    # H^1 = Hom(H_1, Z) + Ext(H_0, Z). H_1=Z, H_0=Z -> Hom(Z,Z) = Z, Ext(Z,Z)=0.
    cohomology[1] = homology[1] # Z
    
    # H^2
    # H^2 = Hom(H_2, Z) + Ext(H_1, Z). H_2=torsion, H_1=Z -> Hom=0, Ext=0.
    cohomology[2] = [] # 0

    # H^i for i>=3
    # H^i = Hom(H_i, Z) + Ext(H_{i-1}, Z). H_i is torsion for i>=2.
    # Hom term is 0. Ext(Torsion, Z) is isomorphic to Torsion.
    # So H^i is isomorphic to H_{i-1}.
    
    # Set a maximum degree for computation based on dimension argument
    max_degree = k
    
    for i in range(3, max_degree + 2):
      h_prev = homology.get(i-1, [])
      cohomology[i] = h_prev
    
    # Find the largest 'a' for which H^a is non-zero
    a = 0
    for i in sorted(cohomology.keys()):
        if cohomology[i]:
            a = i

    result_list = [format_group(cohomology[i]) for i in range(a + 1)]

    print(f"[{', '.join(result_list)}]")


compute_cohomology_list(7)