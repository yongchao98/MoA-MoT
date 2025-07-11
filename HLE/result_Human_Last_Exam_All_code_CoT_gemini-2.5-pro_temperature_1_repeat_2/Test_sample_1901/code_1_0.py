def design_liquid_crystal_molecule():
    """
    This script outlines the design principles for a single-ring liquid crystal
    and generates an example molecular structure based on the provided options.
    """
    print("--- Liquid Crystal Design Report ---")

    # A. Key requirements
    print("\nStep 1: Define Key Requirements")
    print("- Core: Single benzene ring")
    print("- Phase: Nematic or Smectic")
    print("- Transition Temperature: Near room temperature (approx. 20-25 C)")

    # C. Recommended molecular features are the most crucial design step.
    print("\nStep 2: Propose General Molecular Features for Liquid Crystallinity")
    print("- A rigid core is essential for anisotropy. A single benzene ring (Ph) is chosen.")
    print("- To maintain a linear shape, functional groups are attached at para-positions (1 and 4).")
    print("- A flexible terminal group promotes the formation of liquid crystal phases. A flexible alkyl chain (CnH2n+1) is a standard choice.")
    print("- A polar terminal group enhances intermolecular dipole-dipole interactions, which helps stabilize the ordered phases. A nitrile group (-CN) is effective.")

    # E. The general structure
    print("\nStep 3: Construct the General Formula")
    print("Combining these features gives the general molecular structure:")
    print("  AlkylChain - BenzeneRing - PolarGroup")
    print("  C(n)H(2*n+1) - Ph - CN")
    print("Here, 'n' is the number of carbon atoms in the alkyl chain and is the primary variable for tuning properties.")

    # D & F. Specific example and tuning strategy
    print("\nStep 4: Create a Specific Example and Outline Tuning")
    print("To achieve room temperature transitions, a common starting point is a medium-length alkyl chain.")
    # Let's start with n=5, as suggested for the well-known 5CB molecule.
    n = 5
    num_hydrogens = 2 * n + 1
    print(f"We will select a starting value for n = {n}.")
    print(f"Calculating the number of hydrogens: 2 * n + 1 = 2 * {n} + 1 = {num_hydrogens}")

    print("\nFinal Example Molecule:")
    print(f"  - Formula: C{n}H{num_hydrogens}-Ph-CN")
    print(f"  - Name: 1-cyano-4-pentylbenzene")
    print("\nOptimization Strategy:")
    print("  - To adjust the transition temperature, the chain length 'n' can be modified.")
    print("  - Increasing 'n' generally lowers the clearing point (nematic-to-isotropic transition temperature).")
    print("  - Decreasing 'n' generally raises the clearing point.")

    print("\n--- Conclusion ---")
    print("While several options are relevant, Option C provides the most comprehensive and fundamental design guide.")
    print("It explains the 'why' behind the molecular structure by recommending the key features (core, flexible chain, polar group) necessary to induce the desired liquid crystal behavior.")

# Run the design process
design_liquid_crystal_molecule()