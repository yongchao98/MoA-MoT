import sys
import io

# Redirect stdout to capture the output for the final answer format
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_chemistry_problem():
    """
    This function solves the chemical puzzle by calculating molecular formulas and molar masses,
    and by applying chemical principles to deduce the structure of the unknown hydrocarbon X.
    """
    print("Step 1: Determine the molecular formula of substance A1.")
    
    # Elemental composition of A1
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    atomic_masses = {'C': 12.01, 'H': 1.008, 'N': 14.01}

    # Assume 100g of substance to find moles
    moles_c = composition['C'] / atomic_masses['C']
    moles_h = composition['H'] / atomic_masses['H']
    moles_n = composition['N'] / atomic_masses['N']
    
    # Find the smallest mole value to determine the ratio
    min_moles = min(moles_c, moles_h, moles_n)
    
    # Calculate the simplest ratio
    ratio_c = round(moles_c / min_moles)
    ratio_h = round(moles_h / min_moles)
    ratio_n = round(moles_n / min_moles)
    
    print(f"Based on the elemental composition (C: {composition['C']}%, H: {composition['H']}%, N: {composition['N']}%):")
    print(f"The simplest whole number ratio (empirical formula) is approximately C{ratio_c}H{ratio_h}N{ratio_n}.")
    
    # The reaction of a bromo-compound 'A' with excess ammonia to form A1 suggests A1 is a diamine,
    # requiring at least two nitrogen atoms. We test the doubled empirical formula, C4H12N2.
    mf_c, mf_h, mf_n = 4, 12, 2
    molar_mass_A1 = mf_c * atomic_masses['C'] + mf_h * atomic_masses['H'] + mf_n * atomic_masses['N']
    calc_percent_c = (mf_c * atomic_masses['C'] / molar_mass_A1) * 100
    calc_percent_h = (mf_h * atomic_masses['H'] / molar_mass_A1) * 100
    calc_percent_n = (mf_n * atomic_masses['N'] / molar_mass_A1) * 100

    print(f"The molecular formula C{mf_c}H{mf_h}N{mf_n} gives a composition of C: {calc_percent_c:.1f}%, H: {calc_percent_h:.1f}%, N: {calc_percent_n:.1f}%, which matches the given data.")
    print("Therefore, the molecular formula of A1 is C4H12N2.\n")
    
    print("Step 2: Determine the molar mass of the carboxylic acid.")
    mass_acid = 2.16  # g
    vol_koh = 0.030   # L (30 ml)
    conc_koh = 1.0    # M (mol/L)

    # Neutralization is 1:1 for a monoprotic acid
    moles_koh = vol_koh * conc_koh
    moles_acid = moles_koh
    
    molar_mass_acid = mass_acid / moles_acid
    
    print("The neutralization calculation is based on the equation: Molar Mass = mass / (Volume * Concentration).")
    print("Using the provided numbers:")
    print(f"Molar Mass = {mass_acid} g / ({vol_koh} L * {conc_koh} mol/L)")
    print(f"The calculated Molar Mass of the acid is {molar_mass_acid:.0f} g/mol.\n")

    print("Step 3: Deduce the structure of X using the reaction sequence and analytical data.")
    print("The reaction sequence is: X (C4H8 Alkene) -> A (C4H8Br2) -> A1 (C4H12N2) -> A2 (C4H10O2).")
    print("The oxidation of A2 (a C4 diol) produces a carboxylic acid and CO2, indicating the C4 chain is cleaved.")
    
    print("\nWe use the NMR data for A1 ('four types of signals') to identify the correct C4 alkene isomer.")
    print("- If X is but-1-ene, A1 is butane-1,2-diamine. Its four carbons are all chemically different, giving 4 signals. (This is a match)")
    print("- If X is but-2-ene, A1 is butane-2,3-diamine. Due to symmetry, it only has 2 types of carbons. (Not a match)")
    print("- If X is 2-methylpropene, A1 is 2-methylpropane-1,2-diamine. It has 3 types of carbons. (Not a match)\n")
    
    print("The only structure consistent with the NMR data is derived from but-1-ene.")
    
    print("Step 4: Final Conclusion.")
    print("The hydrocarbon X is but-1-ene.")
    print("Full Reaction Pathway:")
    print("X (But-1-ene) + Br2 -> A (1,2-Dibromobutane)")
    print("A + excess NH3 -> A1 (Butane-1,2-diamine)")
    print("A1 + HNO2 -> A2 (Butane-1,2-diol)")
    print("A2 + [O] -> Propanoic acid + CO2")
    
    final_answer = "but-1-ene"
    print(f"\nThe structure of substance X is {final_answer}.")
    
    return final_answer

# Execute the function
solve_chemistry_problem()

# Restore stdout
result_text = captured_output.getvalue()
sys.stdout = old_stdout
print(result_text)

# Extract final answer
final_answer_from_function = result_text.strip().split('\n')[-1].split('is ')[-1].replace('.', '')
print(f"<<<{final_answer_from_function}>>>")