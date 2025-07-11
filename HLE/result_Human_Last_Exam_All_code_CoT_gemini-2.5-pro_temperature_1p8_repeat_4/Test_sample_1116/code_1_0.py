import textwrap

def solve_chemistry_problem():
    """
    This function analyzes the provided chemical reaction and identifies the product.
    """
    # Define the chemical names from the problem
    starting_material = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    reagent = "N-Bromosuccinimide (NBS)"
    product = "2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"

    # Explain the logic step-by-step
    print("### Analysis of the Chemical Reaction ###\n")

    # Step 1: Reactivity Analysis
    print(textwrap.fill("1. The starting material is a large, symmetric molecule with a central electron-accepting core and two pendant thiophene rings. The most reactive positions for electrophilic bromination with NBS are the alpha-positions (C5) of these two outer thiophene rings, as they are the most electron-rich.", 80))
    
    # Step 2: Reaction and Stoichiometry
    print(textwrap.fill("\n2. A complete reaction to brominate both of these equivalent sites requires 2 equivalents of NBS. The use of 2.5 equivalents suggests the reaction might be slow (requiring a slight excess to drive to completion), which is plausible given the deactivating effect of the central acceptor core on the pendant rings.", 80))

    # Step 3: Product Structure
    print(f"\n3. Based on this, the product formed ('the new spot') is the symmetrically di-brominated compound:")
    print(f"\n   Product Name: {product}")
    
    print("\n### Final Equation with Stoichiometry ###\n")
    # As requested, output the numbers in the reaction equation
    print("The balanced chemical equation for the main reaction is:")
    print("1 (Starting Material) + 2 (NBS) ---> 1 (Product) + 2 (Succinimide)")


    # Step 4: Justification with NMR Data
    print("\n### Justification with H-NMR Data ###\n")
    print(textwrap.fill("The key to confirming this structure is the H-NMR data, which shows 'three peaks that are larger than 6.0 ppm'. A simple symmetric product might be expected to show only two peaks. However:", 80))
    print(textwrap.fill("\n- The '2-ethylhexyl' side chains are chiral. This means the product molecule exists as a set of diastereomers (e.g., (R,R), (S,S), and meso-(R,S)).", 80))
    print(textwrap.fill("\n- In the meso-(R,S) diastereomer of the product, the two protons on the central core are no longer chemically equivalent (they are diastereotopic). This results in two separate signals in the NMR spectrum.", 80))
    print(textwrap.fill("\n- The remaining protons on the two outer thiophene rings (at the C3 position) are still equivalent in the meso form and give a single signal.", 80))
    print("\nTherefore, the total number of aromatic signals is 2 (from the core) + 1 (from outer rings) = 3. This matches the experimental data perfectly.")
    
    print("\nConclusion: The new spot is the di-brominated product.")

solve_chemistry_problem()