import dataclasses

@dataclasses.dataclass
class CatalystSystem:
    """Represents a single-site catalyst system."""
    name: str
    metal: str
    ligand: str
    support: str
    polymerization_score: int  # Rated 1-10
    hydrogenolysis_score: int  # Rated 1-10
    rationale: str

def find_optimal_catalyst(candidates):
    """
    Evaluates candidate catalysts to find the one with the best balance
    of polymerization and hydrogenolysis performance.
    A simple scoring model is used that rewards high scores in both categories.
    """
    best_catalyst = None
    max_balanced_score = -1

    print("--- Evaluating Catalyst Candidates ---")
    for cat in candidates:
        # A score that prioritizes systems that are good at both tasks.
        # We multiply the sum by the minimum score to heavily penalize systems
        # that are very poor at one of the two tasks.
        balanced_score = (cat.polymerization_score + cat.hydrogenolysis_score) * min(cat.polymerization_score, cat.hydrogenolysis_score)
        print(f"\nCandidate: {cat.name}")
        print(f"  - Polymerization Score: {cat.polymerization_score}/10")
        print(f"  - Hydrogenolysis Score: {cat.hydrogenolysis_score}/10")
        print(f"  - Calculated Balanced Score: {balanced_score}")

        if balanced_score > max_balanced_score:
            max_balanced_score = balanced_score
            best_catalyst = cat

    return best_catalyst

def main():
    """
    Main function to define candidates and print the optimal result.
    """
    # Define a list of plausible catalyst systems based on scientific literature.
    candidates = [
        CatalystSystem(
            name="Classic Metallocene",
            metal="Zr (Zirconium)",
            ligand="Bis(cyclopentadienyl) 'Cp2'",
            support="Silica (SiO2) + MAO activator",
            polymerization_score=9,
            hydrogenolysis_score=1,
            rationale="Excellent for polymerization but lacks the ability to cleave strong C-C polymer backbones. Not bifunctional."
        ),
        CatalystSystem(
            name="Post-Metallocene (FI Catalyst)",
            metal="Hf (Hafnium)",
            ligand="Bis(phenoxy-imine)",
            support="None (Homogeneous) + Borate activator",
            polymerization_score=10,
            hydrogenolysis_score=2,
            rationale="State-of-the-art for producing specific polymer microstructures. Highly active for polymerization but not designed for degradation."
        ),
        CatalystSystem(
            name="Supported Organometallic Bifunctional Catalyst",
            metal="Zr (Zirconium)",
            ligand="Constrained-Geometry Ligand (e.g., cyclopentadienyl-amido)",
            support="Sulfated Alumina (Al2O3-SO4)",
            polymerization_score=7,
            hydrogenolysis_score=8,
            rationale="A true bifunctional system. The Zr center catalyzes polymerization and hydrogenation. The highly acidic sulfated support facilitates C-C bond cleavage via a carbocation mechanism at high temperatures. The ligand provides stability."
        )
    ]

    optimal_catalyst = find_optimal_catalyst(candidates)

    print("\n--- Optimal Combination Found ---")
    if optimal_catalyst:
        print("The most promising system for both olefin polymerization and polyolefin hydrogenolysis is the Supported Organometallic Bifunctional Catalyst.")
        print("\nFinal Proposed Catalyst Equation:")
        print(f"  1. Metal Center: {optimal_catalyst.metal} (Group 4)")
        print(f"  2. Ligand System: {optimal_catalyst.ligand}")
        print(f"  3. Support Material: {optimal_catalyst.support}")
        print("\nPerformance Profile:")
        print(f"  - Polymerization Efficiency Score: {optimal_catalyst.polymerization_score}")
        print(f"  - Hydrogenolysis Efficiency Score: {optimal_catalyst.hydrogenolysis_score}")
        print("\nJustification:")
        print(optimal_catalyst.rationale)
        print("This combination creates a tandem catalytic effect on a single site, where reaction outcomes can be controlled by tuning temperature and hydrogen pressure.")
    else:
        print("No suitable catalyst could be determined from the provided candidates.")

if __name__ == "__main__":
    main()