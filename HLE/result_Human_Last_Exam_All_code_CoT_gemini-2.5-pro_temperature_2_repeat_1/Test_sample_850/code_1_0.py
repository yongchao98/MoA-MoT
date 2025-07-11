import textwrap

def explain_optical_activity_symmetry():
    """
    Analyzes the symmetry requirements for optical activity in crystal classes
    and addresses the user's specific query.
    """

    # --- Property Definitions based on the 32 Crystallographic Point Groups ---

    # Optical activity is an intrinsic property of chiral materials.
    # Therefore, the set of optically active groups is identical to the set of chiral groups.
    # These groups lack an inversion center, mirror plane, or any roto-inversion axis.
    optically_active_chiral_groups = {
        '1', '2', '3', '4', '6',
        '222', '32', '422', '622',
        '23', '432'
    }

    # All 32 point groups.
    all_32_point_groups = {
        '1', '2', 'm', '-1', '2/m', 'mm2', '222', 'mmm',
        '3', '-3', '32', '3m', '-3m',
        '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm',
        '6', '-6', '6/m', '622', '6mm', '-6m2', '6/mmm',
        '23', 'm-3', '432', '-43m', 'm-3m'
    }

    # Achiral groups are all groups that are not chiral.
    achiral_groups = all_32_point_groups - optically_active_chiral_groups

    # Polar groups possess a unique polar direction.
    polar_groups = {'1', '2', 'm', 'mm2', '3', '3m', '4', '4mm', '6', '6mm'}
    
    # Non-polar groups are all groups that are not polar.
    non_polar_groups = all_32_point_groups - polar_groups

    # --- Analysis of the User's Query ---
    
    print("--- Analysis of the Query: 'achiral', 'non-polar', and 'optically active' ---")
    
    explanation = """
    A fundamental principle in crystallography and optics is that true optical
    activity (the rotation of the plane of polarized light) can only occur in a
    chiral medium. A crystal is chiral if its point group symmetry lacks any
    improper rotations (i.e., no center of inversion, no mirror planes, no
    roto-inversion axes).

    The term 'achiral' means 'not chiral'. An achiral crystal, by definition,
    possesses at least one of these symmetry elements that forbid optical
    activity. Therefore, the conditions 'achiral' and 'optically active'
    are mutually exclusive. A crystal class cannot be both at the same time.
    """
    
    print(textwrap.dedent(explanation).strip())
    
    # We find the intersection of the sets as requested.
    result_set = achiral_groups.intersection(non_polar_groups).intersection(optically_active_chiral_groups)
    
    print("\n[SEARCHING] for crystal classes that are (Achiral AND Non-polar AND Optically Active)...")
    if not result_set:
        print(">>> RESULT: None. The set of such crystal classes is empty, as predicted by physical principles.")
    else:
        # This code block should be unreachable.
        print(f"Result: {sorted(list(result_set))}")

    # --- Providing Related, Valid Information ---
    
    print("\n--- For reference, here are related valid combinations: ---")
    
    # Optically active and non-polar
    optically_active_non_polar = optically_active_chiral_groups.intersection(non_polar_groups)
    print("\n1. Optically Active AND Non-polar crystal classes:")
    print("   (These are chiral and do not have a net dipole)")
    print(f">>> RESULT: {sorted(list(optically_active_non_polar))}")
    
    # Achiral and non-polar
    achiral_non_polar = achiral_groups.intersection(non_polar_groups)
    print("\n2. Achiral AND Non-polar crystal classes:")
    print("   (These are not optically active and do not have a net dipole)")
    print(f">>> RESULT: {sorted(list(achiral_non_polar))}")

if __name__ == "__main__":
    explain_optical_activity_symmetry()

<<<None>>>