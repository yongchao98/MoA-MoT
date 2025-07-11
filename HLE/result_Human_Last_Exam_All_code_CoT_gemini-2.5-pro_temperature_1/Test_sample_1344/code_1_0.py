def propose_dual_function_catalyst():
    """
    This function outlines a proposed optimal combination for a dual-function
    single-site catalyst capable of both olefin polymerization and polyolefin hydrogenolysis.
    The system is based on leading research in plastics upcycling.
    """

    print("--- Proposed Optimal Catalyst System ---")
    print("This system is designed to act as a switchable catalyst:")
    print("- At lower temperatures and without H2, it performs olefin polymerization.")
    print("- At higher temperatures and in the presence of H2, it performs polyolefin hydrogenolysis.\n")

    # The prompt requested outputting numbers in an "equation".
    # We will number the components of the catalyst system.
    # 1. The Metal Center
    # 2. The Ligand System
    # 3. The Support Material

    print("Component 1: The Group IV Metal")
    print("Metal: Zirconium (Zr)")
    print("Reasoning: Zirconium provides an excellent balance of oxophilicity, electrophilicity, and stability. It is a well-established metal for polymerization (e.g., metallocene catalysts) and has shown high activity for C-C bond cleavage when paired with an appropriate ligand and support.\n")

    print("Component 2: The Ligand")
    print("Ligand Class: Constrained-Geometry Ligand (CGC)")
    print("Example: A monocyclopentadienyl-amido ligand, such as [η⁵-C₅Me₄-SiMe₂-N(tBu)]²⁻.")
    print("Reasoning: CGCs create a highly reactive, sterically open, and electronically unsaturated metal center. This geometry is crucial for both allowing monomer access during polymerization and enabling the difficult C-C bond activation steps during hydrogenolysis. The ligand also imparts high thermal stability, which is necessary for the high temperatures required for hydrogenolysis.\n")

    print("Component 3: The Support")
    print("Support: Sulfated Zirconia (SZ)")
    print("Reasoning: This is arguably the most critical component for the dual-functionality. The highly acidic nature of sulfated zirconia is key. It acts as a solid-state activator for the organometallic precursor, generating the cationic active site. Furthermore, the strong Brønsted/Lewis acid sites on the SZ surface are believed to work in tandem with the Zr center to facilitate the cleavage of strong C-C bonds in the polymer chain.\n")

    final_catalyst_description = "A Zirconium-based constrained-geometry single-site catalyst supported on sulfated zirconia (SZ)."
    
    print("--- Final Equation/Combination ---")
    print(f"1. Metal: Zirconium")
    print(f"2. Ligand: Constrained-Geometry Ligand")
    print(f"3. Support: Sulfated Zirconia")
    
    # The final answer in the required format.
    print(f"\n<<<A Zirconium-based constrained-geometry single-site catalyst supported on sulfated zirconia (SZ).>>>")

if __name__ == '__main__':
    propose_dual_function_catalyst()