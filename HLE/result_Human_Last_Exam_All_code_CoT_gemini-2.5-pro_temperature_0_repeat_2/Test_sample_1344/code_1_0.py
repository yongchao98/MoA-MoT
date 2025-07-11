import textwrap

def solve_catalyst_problem():
    """
    Analyzes and determines the optimal combination of a Group IV metal,
    ligand, and support for a dual-function single-site catalyst.
    The catalyst must perform olefin polymerization and polyolefin hydrogenolysis.
    """

    print("### Catalyst Design Analysis ###\n")

    # Step 1: Analyze the dual requirements of the catalyst.
    print("1. The catalyst must perform two key functions:")
    print("   a) Olefin Polymerization: Requires an electrophilic metal center that can coordinate and insert olefin monomers into a growing chain.")
    print("   b) Polyolefin Hydrogenolysis: Requires the ability to cleave strong C-C bonds in the polymer backbone and cap the fragments with hydrogen (H2).\n")

    # Step 2: Select the Group IV Metal.
    print("2. Selecting the Group IV Metal (Ti, Zr, Hf):")
    print("   - Zirconium (Zr) is an excellent choice. It is a well-established metal for high-activity olefin polymerization (e.g., in metallocenes).")
    print("   - Crucially, Zr-hydride species are also known to be highly active for C-C bond activation and hydrogenolysis, making it ideal for the dual role.\n")
    metal = "Zirconium (Zr)"

    # Step 3: Select the Ligand.
    print("3. Selecting the Ligand:")
    print("   - The ligand must provide thermal stability while allowing access for both small olefin monomers and bulky polymer chains.")
    print("   - A 'Constrained-Geometry Catalyst' (CGC) ligand, such as a cyclopentadienyl-amido framework (e.g., [Me2Si(Cp')(NR)]), is optimal.")
    print("   - This ligand's open structure provides high activity and allows the bulky polymer backbone to access the metal center, which is essential for the hydrogenolysis step.\n")
    ligand = "Constrained-Geometry Ligand (e.g., cyclopentadienyl-amido type)"

    # Step 4: Select the Support.
    print("4. Selecting the Support:")
    print("   - A support is needed to create a stable, single-site heterogeneous catalyst.")
    print("   - The support can play a crucial role in C-C bond cleavage. A solid acid support provides a bifunctional mechanism.")
    print("   - Sulfated Zirconia (SZ) is an ideal support. It is a solid superacid that can activate C-C bonds in the polymer backbone, working synergistically with the Zr active site which activates H2 and facilitates the hydrogenation steps.\n")
    support = "Sulfated Zirconia (SZ)"

    # Step 5: Summarize the optimal combination.
    print("### Optimal Catalyst Combination ###\n")
    print(f"Chosen Metal: {metal}")
    print(f"Chosen Ligand Class: {ligand}")
    print(f"Chosen Support: {support}\n")

    final_answer_text = (
        "The optimal combination is a Zirconium (Zr) based Constrained-Geometry Catalyst (CGC) "
        "grafted onto a Sulfated Zirconia (SZ) support. This system is bifunctional: the Zr site "
        "catalyzes polymerization and H2 activation, while the acidic SZ support facilitates the "
        "C-C bond cleavage required for hydrogenolysis."
    )
    
    # Wrapping the text for better display before printing the final answer
    wrapped_text = "\n".join(textwrap.wrap(final_answer_text, width=80))
    
    # Final answer in the required format
    final_answer_formatted = f"<<<{final_answer_text}>>>"
    print(final_answer_formatted)


if __name__ == "__main__":
    solve_catalyst_problem()