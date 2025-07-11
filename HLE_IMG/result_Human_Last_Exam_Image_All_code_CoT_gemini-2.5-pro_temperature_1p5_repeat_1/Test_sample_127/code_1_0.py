import sys

def solve_ring_rotation():
    """
    Analyzes the molecular structure to determine which phenyl rings can rotate freely
    based on steric hindrance rules.
    """
    # Define the structure: pendant rings, their attachment points, and coordination sites.
    # 'metal_bond_pos' for phenylene is the carbon atom number bonded to Iridium.
    # 'metal_bond_pos' for pyridine is 'N'.
    rings = {
        'G': {'type': 'pendant_phenyl', 'attachment_point': {'ring': 'D', 'pos': 6}},
        'H': {'type': 'pendant_phenyl', 'attachment_point': {'ring': 'A', 'pos': 4}},
        'J': {'type': 'pendant_phenyl', 'attachment_point': {'ring': 'C', 'pos': 4}},
        'K': {'type': 'pendant_phenyl', 'attachment_point': {'ring': 'A', 'pos': 6}},
        'L': {'type': 'pendant_phenyl', 'attachment_point': {'ring': 'C', 'pos': 5}},
        'A': {'type': 'phenylene', 'metal_bond_pos': 6},
        'B': {'type': 'pyridine', 'metal_bond_pos': 'N'},
        'C': {'type': 'phenylene', 'metal_bond_pos': 6},
        'D': {'type': 'pyridine', 'metal_bond_pos': 'N'},
    }

    pendant_phenyl_rings = ['G', 'H', 'J', 'K', 'L']
    free_rotators = []

    print("Analyzing steric hindrance for each pendant phenyl ring:")
    print("-" * 50)

    for ring_label in sorted(pendant_phenyl_rings):
        ring_data = rings[ring_label]
        attachment = ring_data['attachment_point']
        parent_ring_label = attachment['ring']
        parent_ring_data = rings[parent_ring_label]
        attachment_pos = attachment['pos']

        is_hindered = False
        reason = "Considered free to rotate as it is not in a sterically hindered position."

        # Rule 1: Hindered if attached to the metallated carbon atom.
        if parent_ring_data.get('type') == 'phenylene' and attachment_pos == parent_ring_data['metal_bond_pos']:
            is_hindered = True
            reason = f"Hindered: Attached directly to the metallated carbon C{attachment_pos} of ring {parent_ring_label}."
        
        # Rule 2: Hindered if attached ortho (adjacent) to the metallated carbon.
        elif parent_ring_data.get('type') == 'phenylene' and abs(attachment_pos - parent_ring_data['metal_bond_pos']) == 1:
            is_hindered = True
            reason = f"Hindered: Attached at C{attachment_pos} of ring {parent_ring_label}, which is ortho to the metallated carbon C{parent_ring_data['metal_bond_pos']}."
        
        # Rule 3: Hindered if attached ortho to the coordinating nitrogen on a pyridine ring.
        elif parent_ring_data.get('type') == 'pyridine' and attachment_pos in [2, 6]:
            is_hindered = True
            reason = f"Hindered: Attached at C{attachment_pos} of pyridine ring {parent_ring_label}, which is ortho to the coordinating Nitrogen."
        
        # Unhindered case (para-position) for clarity
        if not is_hindered and parent_ring_data.get('type') == 'phenylene' and attachment_pos == 4:
            reason = "Free to rotate: Attached at the para-position (C4), which points away from the congested metal center."


        print(f"Ring {ring_label}:")
        print(f"  - Attached to C{attachment_pos} of ring {parent_ring_label}.")
        print(f"  - Verdict: {reason}")
        print("-" * 50)

        if not is_hindered:
            free_rotators.append(ring_label)

    # Sort for consistent output. The problem labels them as J, H in the answer.
    free_rotators.sort(reverse=True)

    print("\nConclusion:")
    print("The phenyl rings that can rotate freely are those with low steric hindrance.")
    # The problem asks to output the final "equation". We will print the letters of the identified rings.
    final_answer_string = ", ".join(free_rotators)
    print(f"Freely rotating rings: {final_answer_string}")

if __name__ == "__main__":
    solve_ring_rotation()