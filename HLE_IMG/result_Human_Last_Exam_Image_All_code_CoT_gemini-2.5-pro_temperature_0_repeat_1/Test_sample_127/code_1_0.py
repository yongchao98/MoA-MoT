import collections

def solve_ring_rotation():
    """
    Analyzes the rotational freedom of substituent phenyl rings in the given complex.
    """
    # Define a named tuple for better readability of ring information
    RingInfo = collections.namedtuple('RingInfo', ['label', 'attachment_pos', 'reason_for_hindrance'])

    # Information about the substituent phenyl rings from the diagram
    substituent_rings = [
        RingInfo(label='G', attachment_pos=6, reason_for_hindrance="ortho to coordinating N, near ligand F"),
        RingInfo(label='H', attachment_pos=4, reason_for_hindrance="para to metalation site, unhindered"),
        RingInfo(label='J', attachment_pos=4, reason_for_hindrance="para to metalation site, unhindered"),
        RingInfo(label='K', attachment_pos=6, reason_for_hindrance="ortho to metalation site, near Ir and ring B"),
        RingInfo(label='L', attachment_pos=6, reason_for_hindrance="ortho to metalation site, near Ir and ring D"),
    ]

    print("Analysis of Phenyl Ring Rotation:")
    print("=" * 60)

    freely_rotating_rings = []

    for ring in substituent_rings:
        # Rule: Rotation is free only if the attachment is at position 4 (para),
        # which is sterically unhindered. Positions 6 (ortho) are hindered.
        is_free = (ring.attachment_pos == 4)
        
        print(f"Ring '{ring.label}':")
        print(f"  - Attachment Position: {ring.attachment_pos} ({'para' if is_free else 'ortho'})")
        print(f"  - Steric Environment: {ring.reason_for_hindrance}")
        print(f"  - Conclusion: Rotation is {'FREE' if is_free else 'RESTRICTED'}")
        print("-" * 40)

        if is_free:
            freely_rotating_rings.append(ring.label)

    # Sort the list for a consistent final answer
    freely_rotating_rings.sort()

    print("\nSummary:")
    print("The rings that can rotate freely are those in sterically unhindered positions.")
    # The prompt asks to "output each number in the final equation".
    # We interpret this as listing the labels of the identified rings.
    print(f"Final set of freely rotating rings: {freely_rotating_rings}")

solve_ring_rotation()