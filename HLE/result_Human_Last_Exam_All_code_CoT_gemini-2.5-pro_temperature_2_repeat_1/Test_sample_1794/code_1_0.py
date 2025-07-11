import math

def solve_phi3_problems():
    """
    This function solves the two user questions regarding 3-loop phi^3 theory.
    """

    # Common parameters
    loop_order = 3
    num_points = 4
    cyclic_ordering = "(1,2,4,3)"

    # --- Part 1: Counting distinct planar graphs ---
    print("--- Part 1: Counting Planar Graphs ---")
    print(f"Theory: Scalar phi^3, L={loop_order} loops, E={num_points} points.")
    print(f"Constraint: Planar for cyclic ordering {cyclic_ordering}, excluding vertex/self-energy corrections (i.e., primitive graphs).")
    
    # In phi^3 theory at 3-loops, 4-points, there are two primitive topologies that can be drawn on a plane:
    # 1. The 'Ladder' graph.
    # 2. The 'Non-Ladder' (or 'Tennis Court') graph.
    # We check which of these is planar for the specified ordering.

    # The Ladder graph is only planar for orderings like (1,2,3,4). It is non-planar for (1,2,4,3).
    is_ladder_planar = False
    
    # The Non-Ladder graph is highly symmetric and is planar for all possible orderings of 4 legs.
    is_nonladder_planar = True
    
    count = 0
    if is_ladder_planar:
        count += 1
    if is_nonladder_planar:
        count += 1
    
    print("\nBased on the topological properties of these graphs:")
    print(f" - Is 'Ladder' graph planar for ordering {cyclic_ordering}? {is_ladder_planar}")
    print(f" - Is 'Non-Ladder' graph planar for ordering {cyclic_ordering}? {is_nonladder_planar}")
    
    final_count = count
    print("\nResult for Part 1:")
    print(f"The number of distinct planar graphs is {final_count}.")

    # --- Part 2: Power of the leading divergence ---
    print("\n--- Part 2: Power of Leading Divergence ---")
    print(f"Process: Massless on-shell 4-point diagram at L={loop_order} loops.")

    print("\nReasoning:")
    print("1. The integrals are Ultraviolet (UV) finite, as the superficial degree of divergence is negative and there are no UV-divergent subgraphs.")
    print("2. However, due to the massless on-shell kinematics, the integrals have Infrared (IR) divergences.")
    print("3. These IR divergences appear as poles in the dimensional regulator epsilon (where d=4-2*epsilon).")
    print("4. For massless theories, the leading IR divergence of an L-loop amplitude typically behaves as (1/epsilon^2)^L.")

    L = loop_order
    # The exponent is 2*L for the term 1/epsilon^(2L). The power is the negative of this.
    power = -2 * L
    
    print("\nCalculation:")
    print(f"The number of loops is L = {L}.")
    print(f"The leading divergence is proportional to 1/epsilon^(2*L).")
    print(f"The power of the leading divergent term is given by the equation: Power = -2 * L")
    print(f"Substituting the value of L: Power = -2 * {L} = {power}")
    
    print("\nResult for Part 2:")
    print(f"The number of the power of the leading divergent term is {power}.")
    
    # This is for the final answer format requested by the user.
    # print(f"\n<<<{final_count}>>>")
    # print(f"<<<{power}>>>")


if __name__ == '__main__':
    solve_phi3_problems()