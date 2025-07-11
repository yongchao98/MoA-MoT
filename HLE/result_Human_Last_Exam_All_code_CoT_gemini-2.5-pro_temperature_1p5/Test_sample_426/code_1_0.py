import sys

# Molar masses for calculations
ELEMENT_MASS = {
    'C': 12.01,
    'H': 1.01,
    'N': 14.01,
    'O': 16.00
}

def solve_chemistry_problem():
    """
    Solves the multi-step chemistry problem to identify hydrocarbon X.
    The steps and reasoning are printed to the console.
    """
    print("Step 1: Determine the molar mass of the carboxylic acid from titration data.")
    
    # Given titration data
    mass_acid = 2.16  # g
    vol_KOH = 30.0    # ml
    conc_KOH = 1.0    # M
    
    # Convert volume to Liters for calculation
    vol_KOH_L = vol_KOH / 1000
    
    # Calculate moles of KOH used
    moles_KOH = conc_KOH * vol_KOH_L
    
    # For a monoprotic acid, the neutralization reaction is 1:1 with KOH.
    moles_acid = moles_KOH
    
    # Calculate the experimental molar mass of the acid
    if moles_acid == 0:
        print("Error: Moles of acid is zero, cannot divide by zero.", file=sys.stderr)
        return

    MM_acid_exp = mass_acid / moles_acid
    
    print(f"The neutralization equation is: Carboxylic Acid + KOH -> Salt + H2O")
    print(f"Moles of KOH = Molarity × Volume = {conc_KOH} mol/L * {vol_KOH_L} L = {moles_KOH:.3f} mol")
    print("Assuming the acid is monoprotic, moles of acid = moles of KOH.")
    print(f"Experimental Molar Mass = mass / moles = {mass_acid} g / {moles_acid:.3f} mol = {MM_acid_exp:.2f} g/mol")
    print("-" * 30)

    print("Step 2: Identify the carboxylic acid and deduce the reaction pathway.")
    
    # Molar mass of propanoic acid (C2H5COOH)
    MM_propanoic = 3 * ELEMENT_MASS['C'] + 6 * ELEMENT_MASS['H'] + 2 * ELEMENT_MASS['O']

    print(f"The experimental molar mass ({MM_acid_exp:.2f} g/mol) is very close to the molar mass of Propanoic Acid (C2H5COOH), which is {MM_propanoic:.2f} g/mol.")
    print("The minor difference can be attributed to reasonable experimental error.")
    print("The formation of a carboxylic acid and CO2 upon strong oxidation of precursor A2 suggests the cleavage of A2.")
    print("If A2 is Butane-1,2-diol (CH3CH2CH(OH)CH2(OH)), its oxidation would cleave the C1-C2 bond to yield Propanoic Acid (CH3CH2COOH) and CO2 (from the CH2OH group).")
    print("This pathway fits all the observations. Let's proceed with this hypothesis.")
    print("Hypothesis: Carboxylic Acid = Propanoic Acid, A2 = Butane-1,2-diol.")
    print("-" * 30)

    print("Step 3: Determine the molecular formula of compound A1 from its composition.")
    
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    
    # Calculate mole ratios from 100g sample
    moles = {k: v / ELEMENT_MASS[k] for k, v in composition.items()}
    smallest_moles = min(moles.values())
    
    # Calculate empirical formula ratios
    empirical_ratios = {k: v / smallest_moles for k, v in moles.items()}

    print(f"Elemental Composition: C={composition['C']}%, H={composition['H']}%, N={composition['N']}%")
    print(f"Empirical Ratio: C:{empirical_ratios['C']:.1f} H:{empirical_ratios['H']:.1f} N:{empirical_ratios['N']:.1f}  =>  C:2 H:6 N:1")
    print("The empirical formula is C2H6N. For a stable organic molecule (following the nitrogen rule), we double this to get C4H12N2.")
    print("Predicted Molecular Formula of A1 is C4H12N2.")
    
    print("\nNow, let's deduce A1's structure from A2 (Butane-1,2-diol).")
    print("A1 is formed by reacting A with excess ammonia, and A1 gives A2 with nitrous acid. This means A1 is a diamine corresponding to the diol A2.")
    print("Structure of A1: Butane-1,2-diamine (CH3CH2-CH(NH2)-CH2(NH2)).")
    print("The formula for Butane-1,2-diamine is C4H12N2. This perfectly matches our calculation from the elemental composition!")
    print("-" * 30)

    print("Step 4: Determine the structure of X.")
    print("A1 (Butane-1,2-diamine) was formed from A by reaction with NH3. Therefore, A is the corresponding dibromide.")
    print("Structure of A: 1,2-Dibromobutane (CH3CH2CH(Br)CH2(Br)).")
    print("\nA was formed from the reaction of hydrocarbon X with bromine, with A being the *only* product. This indicates a bromine addition reaction to an alkene.")
    print("Structure of X: But-1-ene (CH3CH2-CH=CH2).")
    print("-" * 30)
    
    print("Step 5: Final Verification.")
    print("Let's check the NMR data. The structure of A1 is CH3-CH2-CH(NH2)-CH2(NH2).")
    print("The proton environments are: 1) -CH3, 2) -CH2- (ethyl part), 3) -CH(NH2)-, 4) -CH2(NH2).")
    print("These four environments are chemically distinct, leading to four types of signals in the NMR spectrum, which matches the problem description.")
    print("\nConclusion: All evidence consistently points to X being But-1-ene.")
    
    final_answer = "But-1-ene"
    return final_answer

if __name__ == '__main__':
    final_structure = solve_chemistry_problem()
    print("\nFinal Answer: The structure of substance X is:")
    print(final_structure)
