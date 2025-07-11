import pandas as pd

def analyze_mass_spectrum():
    """
    Analyzes the provided mass spectrum data to identify the compound.
    """
    # Step 1: Define the mass spectrum data from the table.
    data = {
        'm/z': [35, 36, 37, 47, 48, 49, 59, 71, 73, 76, 77, 82, 83, 84, 85, 94, 95, 96, 97, 98, 106, 108, 110, 117, 118, 119, 120, 121, 122, 129, 130, 131, 132, 141, 143, 145, 147, 153, 155, 157, 164, 166, 188, 189, 190, 191, 192, 193, 194, 196, 223, 224, 225, 226, 227, 228, 229, 230, 231, 233, 258, 260, 261, 262, 263, 264, 266, 268],
        'I,%': [3.09, 3.38, 1.04, 18.4, 4.62, 6.24, 3.10, 10.6, 3.45, 1.10, 1.08, 4.57, 19.4, 4.01, 7.09, 15.3, 8.79, 10.1, 0.99, 1.04, 10.8, 7.20, 1.25, 2.98, 43.5, 4.41, 29.6, 2.06, 4.88, 2.05, 4.10, 3.42, 1.52, 28.1, 27.8, 9.37, 1.02, 15.9, 16.1, 5.24, 1.36, 1.75, 36.7, 1.61, 48.4, 2.15, 24.5, 1.07, 5.44, 0.45, 60.1, 2.52, 100, 4.40, 66.7, 2.93, 22.2, 0.98, 3.70, 0.95, 17.4, 34.8, 1.53, 29.0, 1.28, 12.9, 3.22, 0.43]
    }
    df = pd.DataFrame(data).set_index('m/z')

    print("Step 1: Identifying the Molecular Ion (M+).")
    # The highest mass cluster of significant intensity starts at m/z 262.
    m_peak = 262
    m_plus_2_peak = 264
    m_plus_4_peak = 266
    
    i_m = df.loc[m_peak]['I,%']
    i_m2 = df.loc[m_plus_2_peak]['I,%']
    i_m4 = df.loc[m_plus_4_peak]['I,%']
    
    print(f"The presumed molecular ion peak (M+) is at m/z = {m_peak} (I={i_m}%).")
    print(f"Isotope peaks are observed at m/z = {m_plus_2_peak} (I={i_m2}%) and m/z = {m_plus_4_peak} (I={i_m4}%).\n")

    print("Step 2: Determining the Presence of Halogens.")
    # The M, M+2, M+4 pattern suggests multiple chlorine atoms.
    # The theoretical intensity ratio for 2 Cl atoms is ~100 : 65 : 10.
    # The observed ratio is 29.0 : 12.9 : 3.22, which normalizes to 100 : 44.5 : 11.1.
    # While not a perfect match, it strongly suggests the presence of two chlorine atoms. Let's verify with fragments.
    print("The pattern of peaks at m/z 262, 264, and 266 suggests the presence of two chlorine atoms.\n")

    print("Step 3: Analyzing Fragmentation.")
    # Find base peak
    base_peak_mz = df['I,%'].idxmax()
    base_peak_intensity = df.loc[base_peak_mz]['I,%']
    print(f"The base peak (most intense) is at m/z = {base_peak_mz} with intensity {base_peak_intensity}%.\n")
    
    # Analyze the [M-Cl]+ fragment
    print("Analysis of the fragment at m/z = 227:")
    loss_from_M = m_peak - 227
    print(f"The peak at m/z=227 corresponds to a loss of {loss_from_M} from the molecular ion (m/z={m_peak}). This matches the mass of a chlorine atom (³⁵Cl).")
    i_227 = df.loc[227]['I,%']
    i_229 = df.loc[229]['I,%']
    ratio_229_227 = i_229 / i_227
    print(f"This fragment has an isotope peak at m/z=229. The intensity ratio I(229)/I(227) is {i_229}/{i_227} = {ratio_229_227:.2f}.")
    print("This ratio of ~0.33 is characteristic of a fragment containing one chlorine atom. This confirms the molecule lost one of its two chlorine atoms to form this fragment.\n")

    print("Analysis of the base peak at m/z = 225:")
    loss_from_227 = 227 - 225
    print(f"The base peak at m/z=225 is 2 amu lighter than the [M-Cl]⁺ fragment at m/z=227.")
    print(f"This corresponds to the loss of a hydrogen molecule (H₂), a common fragmentation process for aromatic/conjugated systems to increase stability.")
    print("Fragmentation pathway: M⁺˙ --(-Cl)--> [M-Cl]⁺ --(-H₂)--> [M-Cl-H₂]⁺")
    print(f"Equation: {m_peak} --(-35)--> {227} --(-2)--> {225}")
    print("This stable [M-Cl-H₂]⁺ ion is the base peak.\n")

    print("Step 4: Deducing the Molecular Formula.")
    # Analyze the [M-2Cl]+ fragment
    mass_of_2cl = 2 * 35
    loss_of_2cl_from_M = m_peak - mass_of_2cl
    print(f"A peak at m/z={loss_of_2cl_from_M} would correspond to the loss of both ³⁵Cl atoms from the molecule.")
    r_core_peak = 192
    i_192 = df.loc[r_core_peak]['I,%']
    print(f"Indeed, there is a peak at m/z = {r_core_peak} (I={i_192}%). This represents the non-halogen core of the molecule, R⁺˙.")
    
    # Deduce formula of R
    mass_of_R = 192
    # C15H12 = 15*12 + 12*1 = 180 + 12 = 192. This is a plausible aromatic/conjugated structure.
    print(f"The mass of this core (R) is {mass_of_R} u. A plausible chemical formula for this is C₁₅H₁₂.")
    print("Therefore, the molecular formula of the original compound is C₁₅H₁₂Cl₂.\n")
    
    print("Step 5: Proposing the Compound and IUPAC Name.")
    print("The fragmentation pattern (loss of Cl, then H₂) suggests a stable, polycyclic aromatic or highly conjugated core.")
    print("A plausible structure for the C₁₅H₁₂ core is methylphenanthrene or phenylindene.")
    print("Assuming a methylphenanthrene core, the compound is a dichloromethylphenanthrene.")
    # We cannot determine the exact isomer positions from the mass spectrum alone.
    # We will provide a name for one possible isomer that fits the derived formula.
    final_name = "9,10-dichloro-1-methylphenanthrene"
    print(f"A possible specific compound is {final_name}.")

    return final_name

final_answer = analyze_mass_spectrum()
print("\n--- Final Answer ---")
print(f"The compound corresponding to the mass spectrum is an isomer of Dichloromethylphenanthrene.")
print(f"A systematic IUPAC name for one such isomer is: {final_answer}")
# Although the exact isomer cannot be determined, this represents the logical conclusion from the spectral data.
# The question asks for *the* compound, which is ambiguous. For the purpose of providing a single answer,
# based on a common C15H12 core, this is a valid candidate.

# After extensive analysis, it seems this is a particularly tricky spectrum, possibly from an advanced source.
# The evidence strongly points to a molecular weight of 262 and a formula of C15H12Cl2.
# The fragmentation pathway is consistent:
# 262 (M+) -> 227 ([M-Cl]+) -> 225 ([M-Cl-H2]+, base peak)
# 262 (M+) -> 192 ([M-2Cl]+)
final_iupac_name = "9,10-dichloro-1-methylphenanthrene"
print(f"\nFinal IUPAC Name: {final_iupac_name}")
print("\nFinal fragmentation equation leading to base peak:")
print("C15H12Cl2+ (m/z=262) -> C15H12Cl+ (m/z=227) -> C15H10Cl+ (m/z=225)")
