import sys

def solve_stereochemistry():
    """
    This script determines the absolute configuration (R/S) of the provided molecule.
    """
    print("Step 1: Identify Stereocenters and IUPAC Numbering")
    print("The molecule is a substituted pentane derivative.")
    print("The principal functional group is the hydroxyl (-OH) group, so it is an alcohol (pentan-2-ol).")
    print("Numbering the longest carbon chain (5 carbons) from left to right to give the -OH group the lowest number (position 2):")
    print("   CH3(1)-CH(OH)(2)-CH(Et)(3)-CH(NH2)(4)-CH3(5)")
    print("The IUPAC name is 4-amino-3-ethylpentan-2-ol.")
    print("There are three stereocenters at positions C2, C3, and C4.")
    print("-" * 30)

    print("Step 2: Determine Configuration of C2")
    print("Groups attached to C2: -OH, -CH(Et)... (C3), -CH3, -H (implied).")
    print("Applying Cahn-Ingold-Prelog (CIP) rules for priority:")
    print("1: -OH (highest atomic number O)")
    print("2: -C3 (carbon attached to two other carbons)")
    print("3: -CH3 (carbon attached to hydrogens)")
    print("4: -H (lowest atomic number)")
    print("3D Orientation: The -OH and -CH3 groups are shown with dashed bonds (away). For a tetrahedral carbon, the implied -H must be on a wedged bond (forward).")
    print("Viewing with the lowest priority group (-H) in front, the path 1 -> 2 -> 3 is clockwise (R).")
    print("Since the lowest priority group points towards us, we REVERSE the configuration.")
    print("Result for C2: R becomes S.")
    print("Configuration at C2 is S.")
    print("-" * 30)

    print("Step 3: Determine Configuration of C3")
    print("Groups attached to C3: -CH(OH)CH3 (C2), -CH(NH2)CH3 (C4), -CH2CH3 (Ethyl), -H (implied).")
    print("Applying CIP rules for priority:")
    print("1: -C2 (C attached to O) > -C4 (C attached to N)")
    print("2: -C4 (C attached to N) > -Ethyl (C attached to C)")
    print("3: -Ethyl")
    print("4: -H")
    print("3D Orientation: The Ethyl group is on a wedged bond (forward). The implied -H must be on a dashed bond (away).")
    print("Viewing with the lowest priority group (-H) in the back, the path 1 -> 2 -> 3 is clockwise.")
    print("Since the lowest priority group points away, the configuration is as it appears.")
    print("Result for C3: Clockwise is R.")
    print("Configuration at C3 is R.")
    print("-" * 30)

    print("Step 4: Determine Configuration of C4")
    print("Groups attached to C4: -NH2, -CH(...) (C3), -CH3, -H (implied).")
    print("Applying CIP rules for priority:")
    print("1: -NH2 (highest atomic number N)")
    print("2: -C3 (C attached to other carbons)")
    print("3: -CH3 (C attached to hydrogens)")
    print("4: -H")
    print("3D Orientation: The -CH3 group is on a dashed bond (away). The main chain C3-C4 bond is in the plane. The C4-NH2 bond is in the plane. The implied -H must be on a wedged bond (forward).")
    print("Viewing with the lowest priority group (-H) in front, the path 1 -> 2 -> 3 is counter-clockwise (S).")
    print("Since the lowest priority group points towards us, we REVERSE the configuration.")
    print("Result for C4: S becomes R.")
    print("Configuration at C4 is R.")
    print("-" * 30)
    
    print("Final Conclusion:")
    print("The absolute configuration of the molecule is determined by combining the configurations of each stereocenter.")
    final_config = "(2S, 3R, 4R)"
    print(f"The absolute configuration is {final_config}.")
    
    # Final answer block
    sys.stdout.write(f"\n<<<{final_config}>>>\n")

solve_stereochemistry()