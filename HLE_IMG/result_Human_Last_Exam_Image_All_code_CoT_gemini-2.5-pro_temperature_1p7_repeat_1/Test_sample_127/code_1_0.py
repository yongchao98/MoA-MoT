def solve_rotation():
    """
    This function determines which phenyl rings in the given molecular structure can rotate freely.
    """
    # Total set of peripheral phenyl rings
    all_peripheral_rings = {'G', 'H', 'J', 'K', 'L'}
    
    # Rings whose rotation is restricted due to significant steric hindrance
    # Rings K and L are at ortho positions (pos 6) to the coordinating C atoms of rings A and C.
    # This proximity to the bulky Iridium coordination sphere prevents free rotation.
    restricted_rings = {'K', 'L'}
    
    # Rings that can rotate freely are the ones not in the restricted set.
    freely_rotating_rings = sorted(list(all_peripheral_rings - restricted_rings))
    
    print("Analysis of Phenyl Ring Rotation:")
    print("-" * 35)
    print("Total peripheral phenyl rings: G, H, J, K, L")
    print("\nStep 1: Identify rings with free rotation.")
    print("Rings G, H, and J are on the periphery, far from the central Iridium atom.")
    print("They are connected by single bonds with low steric hindrance, allowing free rotation.")
    print("\nStep 2: Identify rings with restricted rotation.")
    print("Rings K and L are attached at 'ortho' positions relative to the Ir-C bonds.")
    print("Their rotation is blocked by steric clash with the metal center and other ligands.")
    
    print("\nConclusion:")
    print(f"The rings that can rotate freely are: {', '.join(freely_rotating_rings)}")

solve_rotation()