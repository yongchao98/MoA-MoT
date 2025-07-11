import sys

def solve_feynman_problem():
    """
    This script solves a two-part problem about 3-loop Feynman diagrams
    in massless scalar phi^3 theory.
    """

    # --- Part 1: Number of distinct planar graphs ---
    print("Part 1: Number of distinct planar graphs")
    print("========================================")
    print("The problem asks for the number of distinct planar 3-loop four-point graphs in phi^3 theory with a fixed cyclic ordering of external legs, excluding diagrams with vertex corrections.")
    print("This is standardly interpreted as counting the number of 'primitive' or 'skeleton' graphs, which are graphs that cannot be simplified by removing a self-energy or vertex correction subgraph.")
    
    print("\nFor a 4-point function at 3 loops, it is a known result that there are two such primitive planar topologies:")
    print("1. The Planar Ladder (or Triple-Box) graph.")
    print("2. The Planar Non-Ladder (or Tennis-Court) graph.")
    
    print("\nFor a fixed cyclic ordering of external legs, as given by (1, 2, 4, 3), each of these topologies corresponds to one unique non-isomorphic graph.")
    print("Therefore, the total number of distinct graphs is the sum of the graphs from each topology.")
    
    num_ladder = 1
    num_non_ladder = 1
    total_graphs = num_ladder + num_non_ladder
    
    print(f"\nEquation: {num_ladder} (from ladder topology) + {num_non_ladder} (from non-ladder topology) = {total_graphs}")
    print(f"The number of distinct planar graphs is {total_graphs}.")

    print("\n" + "="*40 + "\n")

    # --- Part 2: Power of the leading divergent term ---
    print("Part 2: Power of the leading divergent term")
    print("============================================")
    print("The problem asks for the power of the leading divergence in the epsilon expansion of the corresponding Feynman integrals near d=4 spacetime dimensions.")
    print("The diagrams are for massless particles with on-shell external momenta (p_i^2 = 0). This setup is known to produce infrared (IR) divergences, which manifest as poles in epsilon (where d = 4 - 2*epsilon).")
    print("The diagrams themselves are primitive and are known to be Ultraviolet (UV) finite. Thus, all divergences are of IR origin.")
    
    print("\nFor massless scattering amplitudes, a general principle states that the leading IR divergence for an L-loop amplitude behaves as 1/epsilon^(2L).")
    print("This leading pole arises from the region of integration where all L loop momenta become 'soft'.")
    
    print("\nIn this problem, the number of loops is L = 3.")
    L = 3
    power_of_divergence = 2 * L
    
    print(f"The calculation for the power of the leading pole is:")
    print(f"Power = 2 * L = 2 * {L} = {power_of_divergence}")
    print(f"The leading divergent term is of the order epsilon^(-{power_of_divergence}).")
    print(f"The number of the power of the leading divergent term is {power_of_divergence}.")

if __name__ == '__main__':
    solve_feynman_problem()