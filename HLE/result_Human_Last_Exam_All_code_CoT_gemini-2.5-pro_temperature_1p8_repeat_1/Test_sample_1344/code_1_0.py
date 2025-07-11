import textwrap

class CatalystSystem:
    """
    A class to represent and describe a proposed dual-function catalyst system.
    """
    def __init__(self, metal, ligand_family, support):
        self.metal = metal
        self.ligand_family = ligand_family
        self.support = support

    def describe(self):
        """Prints a detailed description of the catalyst system and its components."""
        print("### Proposed Optimal Catalyst System for Dual-Functionality ###")
        print("-" * 60)
        
        # --- Metal Component ---
        print(f"Component 1: Metal Center")
        print(f"Choice: {self.metal['name']} ({self.metal['symbol']})")
        reason_metal = textwrap.fill(f"Reasoning: {self.metal['reason']}", width=60)
        print(reason_metal)
        print("-" * 60)

        # --- Ligand Component ---
        print(f"Component 2: Ligand Family")
        print(f"Choice: {self.ligand_family['name']}")
        reason_ligand = textwrap.fill(f"Reasoning: {self.ligand_family['reason']}", width=60)
        print(reason_ligand)
        print("-" * 60)
        
        # --- Support Component ---
        print(f"Component 3: Support Material")
        print(f"Choice: {self.support['name']}")
        reason_support = textwrap.fill(f"Reasoning: {self.support['reason']}", width=60)
        print(reason_support)
        print("-" * 60)

        # --- Overall Rationale ---
        print("Overall Mechanism of Action:")
        overall_rationale = textwrap.fill(
            "This combination creates a 'single-site' catalyst. For "
            "POLYMERIZATION, the active cationic Hf center readily inserts "
            "olefins. For DEGRADATION, adding H₂ and increasing temperature "
            "switches the mechanism. The highly reactive, low-coordinate Hf "
            "center, assisted by the acidic support, can activate the inert C-C "
            "bonds in the polymer chain, leading to hydrogenolysis and cleavage "
            "into shorter, saturated alkanes.",
            width=60
        )
        print(overall_rationale)
        print("-" * 60)

# Define the components of the proposed catalyst system based on literature
optimal_catalyst = CatalystSystem(
    metal={
        'name': "Hafnium",
        'symbol': "Hf",
        'reason': ("Hafnium is a Group IV metal known for forming highly "
                   "active and thermally stable polymerization catalysts. Its complexes "
                   "can be engineered to be sufficiently reactive for C-C bond "
                   "cleavage under hydrogenolysis conditions.")
    },
    ligand_family={
        'name': "Bulky Bis(phenoxy-amine) or Bis(phenoxy-imine) Chelates",
        'reason': ("These 'post-metallocene' ligands are crucial. Their steric bulk "
                   "enforces a specific geometry around the metal center. This "
                   "creates a highly reactive, coordinatively unsaturated site needed for "
                   "cleaving strong C-C bonds, while also providing the stability and "
                   "control required for efficient polymerization.")
    },
    support={
        'name': "Acidic Silica-Alumina (SiO₂-Al₂O₃) or Sulfated Zirconia",
        'reason': ("Immobilizing the Hafnium complex on a solid support is key for "
                   "industrial viability and for breaking down solid plastics. An acidic support "
                   "helps generate and stabilize the highly electrophilic cationic metal "
                   "center required for C-C activation. This creates robust and well-defined single sites.")
    }
)

# Print the final proposed system
optimal_catalyst.describe()