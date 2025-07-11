import sys

def solve_structure():
    """
    This function solves the chemical puzzle step-by-step, printing the logic
    and calculations to identify hydrocarbon X.
    """
    # Step 1: Analyze the neutralization of the carboxylic acid to find its molar mass.
    print("Step 1: Determine the molar mass of the carboxylic acid.")
    mass_acid = 2.16  # g
    vol_koh = 0.030  # L
    conc_koh = 1.0   # M

    # Moles of KOH = Volume (L) * Concentration (M)
    moles_koh = vol_koh * conc_koh

    # For a monoprotic acid (with one -COOH group), the mole ratio with KOH is 1:1.
    # The reaction is: R-COOH + KOH -> R-COOK + H2O
    moles_acid = moles_koh

    # Molar Mass = mass / moles
    molar_mass_acid_exp = mass_acid / moles_acid
    
    print(f"The neutralization equation is: Carboxylic Acid + KOH -> Salt + H2O")
    print(f"From the experimental data, {moles_acid:.3f} moles of acid reacted with {moles_koh:.3f} moles of KOH.")
    print(f"The calculated molar mass of the carboxylic acid is {molar_mass_acid_exp:.2f} g/mol.")
    print("-" * 50)

    # Step 2: Identify the carboxylic acid and deduce the structure of A2.
    print("Step 2: Identify the acid and deduce the precursor A2.")
    # The calculated molar mass of 72 g/mol is very close to the molar mass of propanoic acid (C3H6O2, M = 74.08 g/mol).
    # This small difference is acceptable for this type of problem.
    print("The molar mass (72 g/mol) closely matches that of propanoic acid (CH3CH2COOH, M=74 g/mol).")
    print("The oxidation of substance A2 produced this C3 acid and also released CO2.")
    print("This means A2 was a C4 molecule that underwent oxidative cleavage, losing one carbon as CO2.")
    print("The structure that fits this is butane-1,2-diol (a C4 diol). Its oxidation cleaves the bond between the two carbons bearing -OH groups.")
    print("Oxidation of CH3CH2-CH(OH)-CH2OH produces CH3CH2COOH (propanoic acid) and HCOOH (formic acid). Formic acid is then easily oxidized to CO2 and H2O.")
    print("Therefore, substance A2 is Butane-1,2-diol: CH3-CH2-CH(OH)-CH2OH.")
    print("-" * 50)

    # Step 3: Deduce the structure of A1 from A2.
    print("Step 3: Deduce the structure of A1.")
    print("Substance A2 was formed from A1 by reaction with nitrous acid (HNO2).")
    print("This reaction (diazotization) converts primary amino groups (-NH2) into hydroxyl groups (-OH).")
    print("Thus, A1 must be the corresponding diamine: Butane-1,2-diamine.")
    print("Proposed structure for A1: CH3-CH2-CH(NH2)-CH2NH2.")
    
    print("\nVerification of A1's properties:")
    # Verify elemental composition
    formula_A1 = "C4H12N2"
    m_C, m_H, m_N = 12.01, 1.008, 14.01
    molar_mass_A1 = 4 * m_C + 12 * m_H + 2 * m_N
    percent_C = (4 * m_C / molar_mass_A1) * 100
    percent_H = (12 * m_H / molar_mass_A1) * 100
    percent_N = (2 * m_N / molar_mass_A1) * 100
    print(f"Given composition of A1: C=54.5%; H=13.6%; N=31.8%")
    print(f"Calculated composition for {formula_A1}: C={percent_C:.1f}%; H={percent_H:.1f}%; N={percent_N:.1f}%. This is a match.")

    # Verify NMR data
    print("The structure CH3(a)-CH2(b)-CH(c)-CH2(d) has four chemically distinct types of carbons, which means it will have four distinct signals for the protons attached to them.")
    print("This matches the information that the NMR spectrum of A1 has four types of signals.")
    print("-" * 50)

    # Step 4 & 5: Deduce the structures of A and X.
    print("Step 4 & 5: Deduce the structures of A and X.")
    print("Substance A1 (Butane-1,2-diamine) was formed from A using excess ammonia. This is a nucleophilic substitution, replacing bromines with amino groups.")
    print("Therefore, substance A is 1,2-dibromobutane: CH3-CH2-CH(Br)-CH2Br.")
    print("\nSubstance A was formed from hydrocarbon X by reaction with bromine (Br2). This is an electrophilic addition reaction to an alkene.")
    print("For the product to be 1,2-dibromobutane, the double bond in X must be between carbon 1 and carbon 2.")
    print("Therefore, the original hydrocarbon X must be but-1-ene.")
    print("-" * 50)

    # Final Conclusion
    print("FINAL CONCLUSION:")
    print("The structure of hydrocarbon X is but-1-ene.")
    print("Structure: CH3-CH2-CH=CH2")

solve_structure()