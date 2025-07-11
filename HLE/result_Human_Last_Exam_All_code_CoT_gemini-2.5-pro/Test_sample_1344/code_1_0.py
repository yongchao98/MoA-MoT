import itertools

def generate_catalyst_candidates():
    """
    Generates and prints potential single-site catalyst combinations for
    polyolefin hydrogenolysis and polymerization.

    This script serves as a brainstorming tool for researchers by combining
    common components from the catalysis literature. It does not predict
    the 'optimal' combination but provides plausible starting points for
    further investigation.
    """

    # --- Component Lists ---
    # Group IV metals are standard for this type of catalysis.
    group_iv_metals = ["Zr", "Hf", "Ti"]

    # Ligands are crucial for tuning the catalyst's electronic and steric properties.
    # These represent different classes of ligands used in olefin catalysis.
    ligands = [
        "Cp* (Permethylcyclopentadienyl)",  # Classic, robust metallocene ligand
        "Phenoxy-imine (FI)",              # A well-known post-metallocene ligand class
        "PNP Pincer",                      # Pincer ligands offer high stability
        "N-heterocyclic Carbene (NHC)"     # Strong sigma-donors, create robust catalysts
    ]

    # The support can immobilize the catalyst and influence its activity.
    # 'Homogeneous' means no support is used.
    supports = [
        "Silica (SiO2)",                   # Common, inert support
        "Mesoporous Silica (e.g., SBA-15)",# High surface area support
        "Alumina (Al2O3)",                 # Another common oxide support
        "Homogeneous (No Support)"         # Catalyst is dissolved in the reaction medium
    ]

    print("--- Potential Catalyst System Candidates ---")
    print("This list provides starting points for literature research and experimental design.\n")

    # Generate all unique combinations of the components
    candidates = list(itertools.product(group_iv_metals, ligands, supports))

    for i, (metal, ligand, support) in enumerate(candidates, 1):
        # A representative equation for hydrogenolysis, the process of breaking a
        # carbon-carbon bond with hydrogen. Here, hexane is used as a simple model
        # for a polyethylene chain, being broken down into smaller butane and ethane.
        # The equation is: 1 C6H14 + 1 H2 -> 1 C4H10 + 1 C2H6
        
        print(f"Candidate #{i}:")
        print(f"  Metal:   {metal}")
        print(f"  Ligand:  {ligand}")
        print(f"  Support: {support}")
        
        catalyst_name = f"{metal}/{ligand}/{support}"
        print(f"  Potential Reaction Catalyzed (Hydrogenolysis Model):")
        
        # Print the equation, including the numbers (coefficients) for each molecule.
        print(f"    1 C6H14 + 1 H2 --[{catalyst_name}]--> 1 C4H10 + 1 C2H6")
        
        # Generate a search query to help the user find academic papers.
        search_term = f'"{metal} {ligand} {support} polyolefin hydrogenolysis polymerization"'
        print(f"  Suggested Search Query: {search_term}\n")

if __name__ == '__main__':
    generate_catalyst_candidates()