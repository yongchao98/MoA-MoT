import textwrap

class CatalystComponent:
    """A base class for catalyst components."""
    def __init__(self, name, rationale):
        self.name = name
        self.rationale = rationale

    def describe(self):
        return f"{self.name}\nRationale: {textwrap.fill(self.rationale, 70)}"

def propose_optimal_catalyst():
    """
    This function analyzes chemical principles to propose a hypothetical
    optimal single-site catalyst for both olefin polymerization and
    polyolefin hydrogenolysis.
    """

    # 1. Select the Group IV Metal
    # Hafnium is chosen over Zr and Ti for its superior thermal stability and the
    # typically higher molecular weight polymers it produces, suggesting strong
    # chain propagation ability. Its robust nature is crucial for the harsh
    # conditions often required for C-C bond hydrogenolysis.
    metal = CatalystComponent(
        name="Hafnium (Hf)",
        rationale=("Offers excellent thermal stability, which is essential for "
                   "the high temperatures needed for polymer hydrogenolysis. "
                   "Hf-based catalysts are known for high activity in olefin "
                   "polymerization and can be adapted for bond activation.")
    )

    # 2. Select the Ligand
    # A robust, electron-donating "pincer" ligand is proposed. This architecture
    # provides thermal stability via a tridentate chelate effect. The phosphine
    # arms are strong sigma-donors that can help stabilize the low-valent states
    # potentially needed for H2 activation and hydrogenolysis. The central 'N'
    # can be tuned to modulate the electronic properties and steric bulk.
    ligand = CatalystComponent(
        name="PNP Pincer Ligand (e.g., bis(di-isopropylphosphino)amine)",
        rationale=("The rigid tridentate structure provides high thermal "
                   "stability. The strong sigma-donating phosphine groups can "
                   "support the metal center through multiple oxidation states "
                   "and facilitate challenging bond activations (H-H and C-C), "
                   "while also enabling olefin insertion for polymerization.")
    )

    # 3. Select the Support
    # A Metal-Organic Framework (MOF) provides a highly ordered, high-surface-area
    # support. Using a stable Zr-based MOF like UiO-66 allows for site-isolation of
    # the active Hf catalyst, preventing bimolecular deactivation pathways.
    # The MOF's pores can create a nanoconfined environment that may enhance
    # catalytic selectivity in both reactions.
    support = CatalystComponent(
        name="Zirconium-based Metal-Organic Framework (e.g., UiO-66)",
        rationale=("Provides exceptional thermal and chemical stability. Its high "
                   "surface area and uniform pores allow for the site-isolation of "
                   "the single-site catalyst, preventing deactivation and "
                   "promoting consistent reactivity. The MOF structure itself is "
                   "catalytically inert but provides a crucial stabilizing scaffold.")
    )

    # --- Print the Final Proposed Catalyst System ---
    print("="*60)
    print("PROPOSED BIFUNCTIONAL SINGLE-SITE CATALYST SYSTEM")
    print("="*60)
    print("\nThis system represents a promising research direction for a catalyst capable of\n"
          "both building and breaking down polyolefins.\n")

    print("--- Catalyst Equation ---\n")
    print(f"Component 1 (Metal):   {metal.name}")
    print(f"Component 2 (Ligand):  {ligand.name}")
    print(f"Component 3 (Support): {support.name}")
    print("\n" + "-"*25 + "\n")

    print("--- Component Rationale ---\n")
    print(f"1. Group IV Metal: {metal.describe()}\n")
    print(f"2. Ligand: {ligand.describe()}\n")
    print(f"3. Support: {support.describe()}\n")
    print("="*60)


if __name__ == '__main__':
    propose_optimal_catalyst()