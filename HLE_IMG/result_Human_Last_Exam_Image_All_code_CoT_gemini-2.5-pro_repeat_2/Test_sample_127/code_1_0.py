def solve_molecule_rotation():
    """
    Analyzes the Iridium complex to identify which phenyl rings can rotate freely.
    The primary factor for hindered rotation in this complex is severe steric
    clash, which occurs when a substituent is attached to the same carbon
    atom that is directly bonded (cyclometalated) to the central Iridium ion.
    """

    # Data representation of the phenyl substituents from the diagram.
    # 'parent_ring': The ring the substituent is attached to.
    # 'attachment_carbon': The carbon number on the parent ring for the attachment.
    # 'metalated_carbon': The carbon on the parent ring that is bonded to Iridium.
    rings_info = {
        'G': {'parent_ring': 'C', 'attachment_carbon': 6, 'metalated_carbon': 2},
        'J': {'parent_ring': 'C', 'attachment_carbon': 4, 'metalated_carbon': 2},
        'L': {'parent_ring': 'C', 'attachment_carbon': 5, 'metalated_carbon': 2},
        'K': {'parent_ring': 'A', 'attachment_carbon': 6, 'metalated_carbon': 6},
        'H': {'parent_ring': 'A', 'attachment_carbon': 4, 'metalated_carbon': 6}
    }

    freely_rotating = []
    hindered_ring = ''

    print("Analyzing rotation for each phenyl ring:")
    for ring, info in rings_info.items():
        if info['attachment_carbon'] == info['metalated_carbon']:
            # This is the condition for severe steric hindrance.
            hindered_ring = ring
            print(f" - Ring {ring}: Attached to carbon C{info['attachment_carbon']} of ring {info['parent_ring']}. This carbon is directly bonded to the Iridium. Rotation is HINDERED.")
        else:
            freely_rotating.append(ring)
            print(f" - Ring {ring}: Attached to carbon C{info['attachment_carbon']} of ring {info['parent_ring']}. This is not the metalation site. Rotation is FREE.")

    # Sort for consistent output
    freely_rotating.sort()

    print("\n--- Conclusion ---")
    print(f"The rotation of Ring {hindered_ring} is sterically hindered.")
    print("The rings that can rotate freely are:")
    
    # Print each ring label as requested in the prompt
    for ring in freely_rotating:
        print(ring)

solve_molecule_rotation()
<<<B>>>