import textwrap

def boronic_acid_synthesis_calculator():
    """
    Analyzes the synthesis of (2-bromo-4-chlorophenyl)boronic acid
    and calculates a corrected, more precise reaction stoichiometry.
    """

    # --- Problem Analysis ---
    print("### Problem Analysis ###")
    explanation = """
    The observation of two Boron NMR signals indicates the formation of a side product. The starting material is 2-bromo-4-chloro-1-iodobenzene. The intended reaction is a selective lithium-halogen exchange at the iodine position, followed by quenching with trimethyl borate.

    Desired Reaction:
    Ar-I + n-BuLi -> Ar-Li --(B(OMe)3)--> Ar-B(OH)2
    where Ar = 2-bromo-4-chlorophenyl

    Side Reaction:
    The use of excess n-BuLi (1.05 eq) can cause a second, undesired lithium-halogen exchange at the bromine position, leading to a di-lithiated species. This species forms a di-boronic acid, which is the source of the second NMR signal.

    Ar-I + 2 n-BuLi -> Ar'(Li)2 --(2 B(OMe)3)--> Ar'-(B(OH)2)2
    where Ar' = 4-chlorophenyl

    The solution is to use a very precise, stoichiometric amount (1.00 eq) of n-BuLi to minimize the formation of the di-lithiated side product. The concentration of n-BuLi solutions can decrease over time, so it should be titrated immediately before use for accurate results.
    """
    print(textwrap.dedent(explanation))

    # --- Constants ---
    mw_starting_material = 317.34  # g/mol for 2-bromo-4-chloro-1-iodobenzene
    mw_borate = 103.93             # g/mol for trimethyl borate
    density_borate = 0.932         # g/mL for trimethyl borate

    # --- User Inputs for Calculation ---
    mass_starting_material = 10.0  # g, example scale
    conc_nBuLi = 2.5               # M, typical concentration in hexanes

    # --- Calculations ---
    moles_starting_material = mass_starting_material / mw_starting_material

    # Problematic procedure
    moles_nBuLi_problem = moles_starting_material * 1.05
    vol_nBuLi_problem = (moles_nBuLi_problem / conc_nBuLi) * 1000  # in mL

    # Corrected procedure
    moles_nBuLi_corrected = moles_starting_material * 1.00
    vol_nBuLi_corrected = (moles_nBuLi_corrected / conc_nBuLi) * 1000  # in mL

    # Borate calculation
    moles_borate = moles_starting_material * 5.0
    mass_borate = moles_borate * mw_borate
    vol_borate = mass_borate / density_borate

    # --- Output Final Recommended "Equation" ---
    print("\n### Corrected Stoichiometry Calculation ###\n")
    print(f"Based on {mass_starting_material:.1f} g of starting material and a precise n-BuLi concentration of {conc_nBuLi:.1f} M:\n")

    print("--- Reactants ---")
    print(f"1. 2-bromo-4-chloro-1-iodobenzene (MW={mw_starting_material:.2f} g/mol):")
    print(f"   Mass: {mass_starting_material:.3f} g")
    print(f"   Moles: {moles_starting_material:.4f} mol")

    print(f"\n2. n-BuLi ({conc_nBuLi:.2f} M solution):")
    print(f"   Equivalents: 1.00 eq (Corrected)")
    print(f"   Moles: {moles_nBuLi_corrected:.4f} mol")
    print(f"   Volume: {vol_nBuLi_corrected:.2f} mL")
    print(f"   (Note: The problematic 1.05 eq would have been {vol_nBuLi_problem:.2f} mL. A small measuring error or imprecise concentration can easily lead to excess.)")


    print(f"\n3. Trimethyl Borate (d={density_borate:.3f} g/mL):")
    print(f"   Equivalents: 5.00 eq")
    print(f"   Moles: {moles_borate:.4f} mol")
    print(f"   Mass: {mass_borate:.3f} g")
    print(f"   Volume: {vol_borate:.2f} mL")
    
    print("\n--- Final Equation with Corrected Amounts ---")
    print(f"[{mass_starting_material:.2f} g] 2-bromo-4-chloro-1-iodobenzene + [{vol_nBuLi_corrected:.2f} mL] n-BuLi ({conc_nBuLi:.1f} M) --> Intermediate --(+ [{vol_borate:.2f} mL] B(OMe)3)--> Product")


if __name__ == '__main__':
    boronic_acid_synthesis_calculator()