def solve_feynman_diagram_problems():
    """
    This function provides the answers to the two questions asked by the user.
    
    For the first part, we determine the number of distinct planar graphs
    at 3-loop order for a 4-point scalar interaction, excluding vertex corrections.
    This is a combinatorial problem in quantum field theory. The graphs can be 
    categorized by their structure:

    1. Primitive 3-loop graphs: These cannot be reduced by contracting self-energy
       or vertex subgraphs. There are 2 such topologies.
    2. 1-loop self-energy on the 2-loop primitive graph: This defines 1 unique topology.
    3. 2-loop self-energy on the 1-loop primitive graph: A 2-loop self-energy itself
       has 2 distinct topologies (one primitive, one from nested 1-loop bubbles).
       This gives 2 distinct overall graph topologies.
    4. Two 1-loop self-energies on the 1-loop primitive graph: This configuration
       represents 1 distinct topology.
       
    The total number is the sum of these counts.

    For the second part, we determine the power of the leading divergence term in the
    epsilon expansion of the Feynman integrals for these diagrams. The divergence 
    arises from both Ultraviolet (UV) and Infrared (IR) regions of loop momenta.

    - The maximum UV divergence for an L-loop diagram in this theory is 1/epsilon^L.
      For L=3, this is 1/epsilon^3.
    - The maximum IR divergence for a massless scalar theory is also 1/epsilon^L, as
      only collinear (not soft) divergences are present. For L=3, this is 1/epsilon^3.

    The leading divergence is the highest of these, so the power is 3.
    """

    # Part 1: Number of graphs
    num_primitive_3_loop = 2
    num_se_on_2_loop = 1
    num_2_loop_se_on_1_loop = 2
    num_two_1_loop_se_on_1_loop = 1
    
    total_graphs = (num_primitive_3_loop + 
                    num_se_on_2_loop + 
                    num_2_loop_se_on_1_loop + 
                    num_two_1_loop_se_on_1_loop)
    
    print("--- Part 1: Number of Distinct Planar Graphs ---")
    print(f"Number of primitive 3-loop graphs = {num_primitive_3_loop}")
    print(f"Number of topologies from 1-loop SE on 2-loop primitive = {num_se_on_2_loop}")
    print(f"Number of topologies from 2-loop SE on 1-loop primitive = {num_2_loop_se_on_1_loop}")
    print(f"Number of topologies from two 1-loop SEs on 1-loop primitive = {num_two_1_loop_se_on_1_loop}")
    print(f"Total number of distinct graphs is the sum:")
    print(f"{num_primitive_3_loop} + {num_se_on_2_loop} + {num_2_loop_se_on_1_loop} + {num_two_1_loop_se_on_1_loop} = {total_graphs}")
    print(f"\nFinal answer for Part 1: {total_graphs}")
    
    # Part 2: Power of leading divergence
    leading_divergence_power = 3
    
    print("\n--- Part 2: Power of Leading Divergence ---")
    print("The leading divergence of a 3-loop massless on-shell 4-point diagram in phi^3 theory")
    print("is determined by the most singular behavior in the epsilon expansion (d=4-2*epsilon).")
    print("Both Ultraviolet (UV) and Infrared (IR) divergences contribute.")
    print(f"Maximum UV pole for a 3-loop diagram: 1/epsilon^{leading_divergence_power}")
    print(f"Maximum IR pole for a 3-loop diagram in a scalar theory: 1/epsilon^{leading_divergence_power}")
    print(f"The highest power of the pole 1/epsilon is therefore {leading_divergence_power}.")
    print(f"\nFinal answer for Part 2: {leading_divergence_power}")

solve_feynman_diagram_problems()