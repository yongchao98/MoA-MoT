import math

def is_stable(g, n_special):
    """Checks the stability condition 2g - 2 + n > 0 for a vertex."""
    return 2 * g - 2 + n_special > 0

def calculate_strata():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}.
    This is done by counting stable decorated dual graphs with g=3, n=1, and 2 edges.
    """
    g_total = 3
    n_total = 1
    total_strata = 0

    print("Analyzing codimension 2 boundary strata for M_bar(3,1) by dual graph topology:\n")

    # --- Case A: One vertex, two self-loops ("Figure 8") ---
    # h^1 = 2. Genus condition: g_v + 2 = 3 => g_v = 1.
    # Vertex has 1 leg and 2 loops (each counts as 2 connections).
    # n_v = 1 (leg) + 2 * 2 (loops) = 5.
    count_A = 0
    if is_stable(g=1, n_special=5):
        count_A = 1
    print(f"Case A (1 vertex, 2 loops): Found {count_A} stratum.")
    total_strata += count_A

    # --- Case B: Two vertices, two parallel edges ("Sunset") ---
    # h^1 = 1. Genus condition: g1 + g2 + 1 = 3 => g1 + g2 = 2.
    # Each vertex has 2 edges attached.
    count_B = 0
    # Partition (2,0): Point must be on g=0 component for stability.
    if is_stable(g=2, n_special=2) and is_stable(g=0, n_special=2 + n_total):
        count_B += 1
    # Partition (1,1): Point on one component (the other is then distinct but stable).
    if is_stable(g=1, n_special=2 + n_total) and is_stable(g=1, n_special=2):
        count_B += 1
    print(f"Case B (2 vertices, 2 parallel edges): Found {count_B} strata.")
    total_strata += count_B

    # --- Case C: Three vertices in a chain ---
    # h^1 = 0. Genus condition: g1 + g2 + g3 = 3.
    # Terminal vertices (v1, v3) must have g > 0 to be stable, as they have only 1 edge attached.
    # This forces the genus partition to be (1,1,1).
    count_C = 0
    # Point on an end vertex (v1 or v3 are symmetric positions).
    # v1(marked): g=1, n=1+1=2. v2: g=1, n=2. v3: g=1, n=1.
    if is_stable(1, 2) and is_stable(1, 2) and is_stable(1, 1):
        count_C += 1
    # Point on the middle vertex.
    # v1: g=1, n=1. v2(marked): g=1, n=2+1=3. v3: g=1, n=1.
    if is_stable(1, 1) and is_stable(1, 3) and is_stable(1, 1):
        count_C += 1
    print(f"Case C (3-vertex chain): Found {count_C} strata.")
    total_strata += count_C
    
    # --- Case D: Two vertices, one with a loop and one as a tail ---
    # h^1 = 1. Genus condition: g1 + g2 + 1 = 3 => g1 + g2 = 2.
    # The tail vertex (v2) must have g > 0 for stability.
    # Possible partitions (g_loop, g_tail): (1,1) and (0,2).
    count_D = 0
    # Partition (g_loop=1, g_tail=1)
    # Point on loop vertex: v_loop(g=1, n=3+1=4), v_tail(g=1, n=1)
    if is_stable(1, 4) and is_stable(1, 1):
        count_D += 1
    # Point on tail vertex: v_loop(g=1, n=3), v_tail(g=1, n=1+1=2)
    if is_stable(1, 3) and is_stable(1, 2):
        count_D += 1
    # Partition (g_loop=0, g_tail=2)
    # Point on loop vertex: v_loop(g=0, n=3+1=4), v_tail(g=2, n=1)
    if is_stable(0, 4) and is_stable(2, 1):
        count_D += 1
    # Point on tail vertex: v_loop(g=0, n=3), v_tail(g=2, n=1+1=2)
    if is_stable(0, 3) and is_stable(2, 2):
        count_D += 1
    print(f"Case D (loop-and-tail): Found {count_D} strata.")
    total_strata += count_D
    
    print("\n--- Summary ---")
    print(f"The total number of codimension 2 boundary strata is the sum of counts from each case:")
    print(f"{count_A} + {count_B} + {count_C} + {count_D} = {total_strata}")

calculate_strata()