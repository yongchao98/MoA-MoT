import textwrap

def suggest_catalyst_system():
    """
    Identifies and describes an optimal Group IV single-site catalyst system
    for both olefin polymerization and polyolefin hydrogenolysis.
    """

    # 1. Define the optimal components based on chemical principles
    metal = {
        "group": "IV",
        "symbol": "Zr",
        "name": "Zirconium",
        "reason": "Excellent balance of high activity for both polymerization and C-C bond cleavage."
    }

    ligand = {
        "name": "bis(pentamethylcyclopentadienyl)",
        "abbreviation": "Cp*2",
        "formula_fragment": "(C5(CH3)5)2",
        "reason": "Provides exceptional thermal stability required for high-temperature hydrogenolysis."
    }

    precatalyst_form = {
        "name": "Dimethylbis(pentamethylcyclopentadienyl)zirconium",
        "formula": "Cp*2ZrMe2",
        "reason": "A stable precatalyst from which the active site is generated."
    }

    support = {
        "name": "Sulfated Zirconia",
        "formula": "ZrO2/SO4^2-",
        "reason": ("Acts as a solid superacid activator. It abstracts a methyl group to create the "
                   "highly electrophilic active site, [Cp*2ZrMe]+, and its acidity assists in "
                   "C-C bond cleavage at high temperatures.")
    }

    # 2. Print the detailed description of the final system
    print("--- Optimal Single-Site Catalyst System ---")
    print("\nBased on current research for a dual-function catalyst, the suggested optimal combination is:")
    print(f"\nCATALYST PRECURSOR: {precatalyst_form['name']} ({precatalyst_form['formula']})")
    print(f"SUPPORT / ACTIVATOR: {support['name']} ({support['formula']})")

    print("\n--- Component Breakdown ---")

    # Metal Details
    print(f"\n1. METAL: {metal['name']} ({metal['symbol']})")
    print(textwrap.fill(f"   - Rationale: {metal['reason']}", width=80))

    # Ligand Details
    print(f"\n2. LIGAND: {ligand['name']} ({ligand['abbreviation']})")
    print(textwrap.fill(f"   - Rationale: {ligand['reason']}", width=80))

    # Support Details
    print(f"\n3. SUPPORT/ACTIVATOR: {support['name']}")
    print(textwrap.fill(f"   - Rationale: {support['reason']}", width=80))

    # 4. Explain the dual functionality
    print("\n--- Dual-Function Mechanism ---")
    print(textwrap.fill("A. For Olefin Polymerization (Lower Temp): The supported catalyst forms the "
                        "[Cp*2Zr-R]+ active site, which readily polymerizes olefins to high molecular "
                        "weight polyolefins.", width=80))
    print(textwrap.fill("\nB. For Polyolefin Hydrogenolysis (Higher Temp, H2 pressure): The same thermally "
                        "stable active site, in concert with the acidic support, facilitates the cleavage "
                        "of C-C bonds in the polymer backbone, followed by hydrogenation to produce short "
                        "saturated alkanes (e.g., diesel and gasoline range hydrocarbons).", width=80))

if __name__ == '__main__':
    suggest_catalyst_system()