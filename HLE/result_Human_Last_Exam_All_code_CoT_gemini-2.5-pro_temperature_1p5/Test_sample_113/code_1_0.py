def count_codimension_2_strata():
    """
    Calculates the number of codimension 2 boundary strata of M_3,1.
    This corresponds to counting stable dual graphs with total genus 3,
    1 marked point, and 2 edges.
    """
    total_genus = 3
    num_marked_points = 1
    
    # Store the counts for each type of graph
    counts = {}

    print("Analyzing codimension 2 boundary strata of M_{3,1}...\n")
    print("This is equivalent to counting stable dual graphs with total genus 3, 1 marked point, and 2 edges.\n")
    
    # --- Type 1: One vertex with two loops ---
    # This represents an irreducible curve of genus 3 with two nodes.
    # The vertex has genus 3, one marked point, and two loops.
    # Each loop corresponds to 2 "attachments" or special points.
    print("--- Type 1: Irreducible curve with 2 nodes (1 vertex, 2 loops) ---")
    g = 3
    # Number of special points (k): 1 marked point + 2 nodes (each node adds 1 to k)
    # Alternatively, using dual graph: 1 leg + 2 loops * 2 half-edges = 5 attachments
    # Let's use k = number of nodes + number of marked points for clarity.
    k = 2 + num_marked_points
    if 2 * g - 2 + k > 0:
        counts['type1'] = 1
        print(f"Configuration: Genus=3, 1 m.p. is stable (2*3-2+3 = 7 > 0).")
    else:
        counts['type1'] = 0
    print(f"Number of Type 1 strata: {counts['type1']}\n")


    # --- Type 2: Two vertices, two parallel edges ---
    print("--- Type 2: Two components joined at 2 nodes (2 vertices, 2 parallel edges) ---")
    type2_strata = []
    genus_partitions = [(0, 3), (1, 2)] # g1+g2=3, g1<=g2 to avoid duplicates
    for g1, g2 in genus_partitions:
        # Case A: marked point on the first component (g1)
        k1 = 2 + num_marked_points
        k2 = 2
        if 2 * g1 - 2 + k1 > 0 and 2 * g2 - 2 + k2 > 0:
            desc = f"(g1,g2)=({g1},{g2}), m.p. on g1 component"
            type2_strata.append(desc)
            print(f"Found stable stratum: {desc}")
        
        # Case B: marked point on the second component (g2)
        k1 = 2
        k2 = 2 + num_marked_points
        if 2 * g1 - 2 + k1 > 0 and 2 * g2 - 2 + k2 > 0:
            desc = f"(g1,g2)=({g1},{g2}), m.p. on g2 component"
            type2_strata.append(desc)
            print(f"Found stable stratum: {desc}")
    counts['type2'] = len(type2_strata)
    print(f"Number of Type 2 strata: {counts['type2']}\n")
    
    
    # --- Type 3: Two vertices, one edge, one loop ---
    print("--- Type 3: A component with a node connected to another self-nodal component ---")
    # v_e: end vertex (1 edge), v_l: loop vertex (1 edge, 1 loop)
    type3_strata = []
    genus_permutations = [(0, 3), (3, 0), (1, 2), (2, 1)] # All (g_e, g_l) where g_e+g_l=3
    for g_e, g_l in genus_permutations:
        # Case A: marked point on end vertex
        k_e = 1 + num_marked_points; k_l = 2
        if 2*g_e-2+k_e > 0 and 2*g_l-2+k_l > 0:
            desc = f"End genus {g_e} (with m.p.), Loop genus {g_l}"
            type3_strata.append(desc)
            print(f"Found stable stratum: {desc}")
            
        # Case B: marked point on loop vertex
        k_e = 1; k_l = 2 + num_marked_points
        if 2*g_e-2+k_e > 0 and 2*g_l-2+k_l > 0:
            desc = f"End genus {g_e}, Loop genus {g_l} (with m.p.)"
            type3_strata.append(desc)
            print(f"Found stable stratum: {desc}")
    counts['type3'] = len(type3_strata)
    print(f"Number of Type 3 strata: {counts['type3']}\n")


    # --- Type 4: Three vertices in a chain ---
    print("--- Type 4: Three components in a chain ---")
    type4_count = 0
    
    # Case A: marked point on an end component (e.g., v1).
    # Symmetry means m.p. on v1 is same type as m.p. on v3.
    # Stability requires g1,g2,g3 >= 1. The only partition is (1,1,1).
    type4_count += 1
    print("Found stable stratum: Genera (1,1,1), m.p. on an end component.")
    
    # Case B: marked point on the middle component (v2).
    # Stability requires g1,g3 >= 1 and g2 >= 0.
    middle_configs = [(1, 1, 1), (1, 0, 2), (2, 0, 1)]
    for g1, g2, g3 in middle_configs:
         type4_count += 1
         print(f"Found stable stratum: Genera ({g1},{g2},{g3}), m.p. on the middle component.")
    counts['type4'] = type4_count
    print(f"Number of Type 4 strata: {counts['type4']}\n")

    # --- Final Summation ---
    print("--- Total Count ---")
    total = sum(counts.values())
    calculation = " + ".join(str(c) for c in counts.values())
    print(f"The total number of strata is the sum of the counts from each type:")
    print(f"Total = {calculation} = {total}")
    
    return total

if __name__ == "__main__":
    count_codimension_2_strata()