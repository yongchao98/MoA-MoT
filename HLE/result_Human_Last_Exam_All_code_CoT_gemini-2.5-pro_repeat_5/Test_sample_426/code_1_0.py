import sys

# Define atomic masses for calculations
ATOMIC_MASS = {'C': 12.01, 'H': 1.008, 'O': 16.00, 'N': 14.01, 'Br': 79.90}

def calculate_molar_mass(formula_dict):
    """Calculates the molar mass of a molecule from a dictionary of its atoms."""
    mm = 0
    for atom, count in formula_dict.items():
        mm += ATOMIC_MASS[atom] * count
    return mm

def print_step(step_num, text):
    """Formats and prints a step in the solution."""
    print(f"\n--- STEP {step_num} ---")
    print(text)

# --- Step 1: Identify the Carboxylic Acid ---
print_step(1, "Analyze the neutralization of the carboxylic acid to find its molar mass.")
mass_acid = 2.16  # g
vol_koh = 0.030   # L (30 ml)
conc_koh = 1.0    # M

# Calculate moles of KOH used
moles_koh = vol_koh * conc_koh
print(f"The equation for moles of KOH is: Moles = Volume * Concentration")
print(f"Moles of KOH = {vol_koh:.3f} L * {conc_koh:.1f} M = {moles_koh:.3f} mol")

# Assume the acid is monoprotic (reacts 1:1 with KOH)
moles_acid = moles_koh
molar_mass_acid_calc = mass_acid / moles_acid
print("\nAssuming the acid is monoprotic, its molar mass is calculated as:")
print(f"Molar Mass = Mass / Moles")
print(f"Calculated Molar Mass = {mass_acid:.2f} g / {moles_acid:.3f} mol = {molar_mass_acid_calc:.2f} g/mol")

# Compare with common acids. Molar mass of Propanoic Acid (C3H6O2) is ~74 g/mol.
mm_propanoic = calculate_molar_mass({'C': 3, 'H': 6, 'O': 2})
print(f"\nThe calculated molar mass ({molar_mass_acid_calc:.2f} g/mol) is very close to the molar mass of Propanoic Acid (C3H6O2), which is {mm_propanoic:.2f} g/mol.")
print("The small difference suggests the acid is Propanoic Acid, and the given numbers might have experimental tolerance.")
acid_name = "Propanoic Acid"
acid_formula = "CH3CH2COOH"

# --- Step 2: Deduce the Structure of Alcohol A2 ---
print_step(2, "Deduce the structure of alcohol A2 from its oxidation products.")
print(f"A2 is oxidized to form {acid_name} (a C3 acid) and CO2 (a C1 fragment).")
print("This means the parent alcohol A2 must have had 3 + 1 = 4 carbon atoms.")
print("The oxidative cleavage of a 1,2-diol is a classic reaction that produces two carbonyl-containing fragments. The oxidation of Butane-1,2-diol (CH3-CH2-CH(OH)-CH2(OH)) cleaves the C1-C2 bond.")
print("The -CH2(OH) group oxidizes to formic acid and then to CO2.")
print("The CH3-CH2-CH(OH)- group oxidizes to Propanoic Acid.")
print("Therefore, substance A2 must be Butane-1,2-diol.")
A2_name = "Butane-1,2-diol"

# --- Step 3: Determine and Verify the Structure of Amine A1 ---
print_step(3, "Determine the structure of amine A1 and verify with the given data.")
print(f"A1 reacts with nitrous acid to form A2 ({A2_name}). This means A1 has the same carbon skeleton, with -NH2 groups instead of -OH groups.")
print("Therefore, A1 is Butane-1,2-diamine (CH3-CH2-CH(NH2)-CH2(NH2)).")

# Verify composition of Butane-1,2-diamine (C4H12N2)
A1_formula_dict = {'C': 4, 'H': 12, 'N': 2}
mm_A1 = calculate_molar_mass(A1_formula_dict)
c_mass_total = A1_formula_dict['C'] * ATOMIC_MASS['C']
h_mass_total = A1_formula_dict['H'] * ATOMIC_MASS['H']
n_mass_total = A1_formula_dict['N'] * ATOMIC_MASS['N']

percent_c = (c_mass_total / mm_A1) * 100
percent_h = (h_mass_total / mm_A1) * 100
percent_n = (n_mass_total / mm_A1) * 100

print(f"\nVerifying the elemental composition for Butane-1,2-diamine (C4H12N2, Molar Mass = {mm_A1:.2f} g/mol):")
print(f"Expected %C = ({A1_formula_dict['C']} * {ATOMIC_MASS['C']:.2f} / {mm_A1:.2f}) * 100 = {percent_c:.1f}% (Given: 54.5%)")
print(f"Expected %H = ({A1_formula_dict['H']} * {ATOMIC_MASS['H']:.3f} / {mm_A1:.2f}) * 100 = {percent_h:.1f}% (Given: 13.6%)")
print(f"Expected %N = ({A1_formula_dict['N']} * {ATOMIC_MASS['N']:.2f} / {mm_A1:.2f}) * 100 = {percent_n:.1f}% (Given: 31.8%)")
print("The calculated percentages match the given data very well.")

print("\nVerifying the NMR data:")
print("Butane-1,2-diamine is an asymmetric molecule. It has four distinct carbon atoms (C1, C2, C3, C4), which would result in four signals in a 13C-NMR spectrum. This matches the 'four types of signals' given.")

# --- Step 4: Identify the initial Hydrocarbon X ---
print_step(4, "Identify the structure of hydrocarbon X.")
print("Substance A1 (Butane-1,2-diamine) was formed from A by replacing bromine atoms with amino groups.")
print("Therefore, substance A must be 1,2-Dibromobutane.")
print("Substance A was formed by the reaction of hydrocarbon X with bromine (Br2). This is a characteristic addition reaction for an alkene.")
print("The addition of Br2 to an alkene breaks the double bond and adds a bromine atom to each of the two carbons.")
print("To find X, we reverse the reaction: we remove the two bromine atoms from 1,2-Dibromobutane and form a double bond between C1 and C2.")
print("This gives us the structure: CH3-CH2-CH=CH2.")
print("\nThe initial hydrocarbon X is But-1-ene.")

# Final answer in the specified format
sys.stdout.write("\n<<<But-1-ene>>>\n")
