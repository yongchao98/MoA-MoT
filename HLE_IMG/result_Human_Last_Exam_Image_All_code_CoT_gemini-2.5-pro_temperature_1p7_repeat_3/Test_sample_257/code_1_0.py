import textwrap

def solve_nmr_problem():
    """
    Analyzes the reaction and determines the 1H NMR splitting pattern and
    integration for the most deshielded proton of Compound 1.
    """
    analysis_steps = [
        "### Step-by-step Reasoning:",
        "",
        "1. **Analysis of the Starting Material (Pr-DAOTA):**",
        "   - The starting molecule is a cationic polycyclic aromatic compound. In such systems, protons attached to the aromatic rings are significantly deshielded (appear at high ppm values).",
        "   - The most deshielded proton is typically the isolated one on the central aromatic ring, located between two fused ring systems. This proton is in a unique electronic environment, being part of a positively charged core.",
        "   - Looking at its local environment, this central proton has no protons on adjacent carbons (no vicinal neighbors). According to the n+1 rule, its signal will not be split.",
        "   - **Conclusion for Starting Material:** The most deshielded peak is a singlet (s) with a relative integration of 1H.",
        "",
        "2. **Analysis of the Reaction:**",
        "   - The reaction uses concentrated sulfuric acid (conc. H₂SO₄), a strong acid and sulfonating agent.",
        "   - The product, Compound 1, is described as 'water-soluble'. The introduction of a sulfonic acid group (-SO₃H) onto an organic molecule is a standard method to increase its water solubility.",
        "   - Therefore, the reaction is an electrophilic aromatic sulfonation, where one or more -SO₃H groups are added to the aromatic rings of Pr-DAOTA.",
        "",
        "3. **Analysis of the Product (Compound 1):**",
        "   - Compound 1 is the sulfonated version of Pr-DAOTA. The sulfonation occurs on the outer aromatic rings, not the central one.",
        "   - The unique central proton is still present in the structure of Compound 1. The sulfonation reaction does not affect its position or its neighbors.",
        "   - As this proton still has no protons on adjacent carbons, its splitting pattern remains a singlet.",
        "   - Since it's still a single, unique proton, its integration remains 1H.",
        "   - This proton is expected to remain the most deshielded in the molecule, as the factors making it highly deshielded are inherent to the core cationic structure.",
        ""
    ]

    for step in analysis_steps:
        print(textwrap.fill(step, width=90))

    splitting_pattern = "singlet"
    integration = "1H"

    print("### Final Answer:")
    print(f"The splitting pattern of the highest deshielded proton peak in Compound 1 is: {splitting_pattern}")
    print(f"The integration of this peak is: {integration}")

solve_nmr_problem()