import math

def calculate_phi4_symmetry_factors():
    """
    Calculates the symmetry factors for each second-order vacuum bubble diagram
    in phi^4 scalar field theory and sums them.
    """

    # Diagram 1: Figure-eight diagram
    # This diagram has two vertices connected by four propagators.
    # Its symmetry factor is the product of the symmetry from permuting vertices
    # and the symmetry from permuting propagators.
    
    # Symmetry from swapping the 2 identical vertices
    vert_sym_d1 = 2
    # Symmetry from permuting the 4 identical propagators between the vertices
    prop_sym_d1 = math.factorial(4)
    s1 = vert_sym_d1 * prop_sym_d1

    # Diagram 2: "Saturn" diagram
    # This diagram has two vertices connected by two propagators, and each
    # vertex has one self-loop.
    
    # Symmetry from swapping the 2 identical vertices
    vert_sym_d2 = 2
    # Symmetry from permuting the 2 identical propagators between the vertices
    prop_sym_d2 = math.factorial(2)
    # Symmetry from the two self-loops (each loop contributes a factor of 2)
    loop_sym_d2 = 2 * 2
    s2 = vert_sym_d2 * prop_sym_d2 * loop_sym_d2

    # Summing the symmetry factors
    total_s = s1 + s2

    print("Calculation of symmetry factors for second-order vacuum bubbles in phi^4 theory:")
    print("-" * 70)
    print("1. Figure-eight diagram:")
    print(f"   Symmetry factor S1 = (vertex permutations) * (propagator permutations)")
    print(f"   S1 = {vert_sym_d1} * {prop_sym_d1} = {s1}")
    print("\n2. 'Saturn' diagram:")
    print(f"   Symmetry factor S2 = (vertex perm.) * (propagator perm.) * (loop factors)")
    print(f"   S2 = {vert_sym_d2} * {prop_sym_d2} * {loop_sym_d2} = {s2}")
    print("-" * 70)
    print("Total sum of symmetry factors:")
    print(f"{s1} + {s2} = {total_s}")

calculate_phi4_symmetry_factors()