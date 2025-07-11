import sys
import io

# Redirect stdout to capture the print output for the final answer format
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_compound_and_determine_iupac_name():
    """
    Analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    The analysis is printed step-by-step.
    """

    print("Step 1: Analysis of Mass Spectrum (MS)")
    mw = 135
    print(f"The molecular ion peak (M+) is observed at m/z {mw}.")
    print(f"The odd molecular weight of {mw} indicates the presence of an odd number of nitrogen atoms (Nitrogen Rule), most likely one.")
    print("A significant fragment peak at m/z 91 is characteristic of a benzyl cation ([C6H5CH2]+), suggesting a benzyl or substituted benzyl group.")
    print("If the formula is C9H13N, the calculated molecular weight is (9 * 12.01) + (13 * 1.008) + (1 * 14.007) = 135.21. This matches the M+ peak.")
    print("Therefore, the molecular formula is likely C9H13N.")
    print("-" * 50)

    print("Step 2: Degree of Unsaturation (DoU) Calculation")
    C, H, N = 9, 13, 1
    # DoU = C - H/2 + N/2 + 1
    dou = C - (H/2) + (N/2) + 1
    print(f"For the formula C{C}H{H}N{N}, the calculation is: DoU = {C} - ({H}/2) + ({N}/2) + 1 = {int(dou)}")
    print(f"A Degree of Unsaturation of {int(dou)} strongly suggests the presence of a benzene ring (which accounts for 4 DoU).")
    print("-" * 50)

    print("Step 3: Analysis of Infrared (IR) Spectrum")
    print("- Two absorption bands in the 3300-3400 cm-1 region indicate a primary amine (-NH2) group.")
    print("- Absorption bands above 3000 cm-1 indicate aromatic C-H stretching.")
    print("- Absorption bands below 3000 cm-1 indicate aliphatic C-H stretching.")
    print("- Absorption bands at ~1600 and ~1450 cm-1 indicate C=C stretching within an aromatic ring.")
    print("Conclusion from IR: The molecule contains a primary amine and a benzene ring.")
    print("-" * 50)

    print("Step 4: Analysis of 13C NMR and DEPT-135 Spectra")
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    print(f"The 13C NMR spectrum shows 7 signals at shifts (ppm): {c13_shifts}.")
    print("DEPT-135 information indicates one negative signal (one CH2 group) and five positive signals (CH and CH3 groups).")
    print("Signal Assignment:")
    print(f"  - 145.1 ppm: Aromatic quaternary carbon (Cq), which does not appear in DEPT-135.")
    print(f"  - 128.5, 127.3, 126.3 ppm: Three signals for the aromatic CH carbons.")
    print(f"  - 43.5 ppm: Aliphatic CH2 group (negative DEPT-135 signal). Likely the benzylic -CH2-.")
    print(f"  - 49.6 ppm: Aliphatic CH group (positive DEPT-135 signal). Likely the carbon bearing the -NH2 group.")
    print(f"  - 19.2 ppm: Aliphatic CH3 group (positive DEPT-135 signal).")
    print("Summary of carbons: 1 Cq, ~5 CH (aromatic), 1 CH2, 1 CH (aliphatic), 1 CH3. Total = 9 carbons. This is consistent with the formula C9H13N.")
    print("-" * 50)

    print("Step 5: Analysis of 1H NMR and HSQC Spectra")
    print("- 1H NMR ~7.2 ppm (multiplet, 5H): Protons on a monosubstituted benzene ring (C6H5-).")
    print("- 1H NMR ~1.1 ppm (doublet, 3H): A methyl group (-CH3) adjacent to a CH group.")
    print("- 1H NMR ~2.7-3.1 ppm (complex multiplets, 3H total): Corresponds to a CH2 group and a CH group.")
    print("- The HSQC spectrum confirms these C-H connections:")
    print("  - C at 19.2 ppm correlates with H at ~1.1 ppm (-CH3 group).")
    print("  - C at 43.5 ppm correlates with H at ~2.7 ppm (-CH2- group).")
    print("  - C at 49.6 ppm correlates with H at ~3.1 ppm (-CH- group).")
    print("-" * 50)

    print("Step 6: Structure Assembly and IUPAC Naming")
    print("All data collectively point to a single structure. The fragments are a phenyl group, a CH2, a CH, a CH3, and an NH2 group.")
    print("The connectivity is determined as C6H5-CH2-CH(NH2)-CH3.")
    print("\nTo determine the IUPAC name:")
    parent_chain = 3
    amine_position = 2
    phenyl_position = 1
    print(f"1. The longest carbon chain containing the principal functional group (amine) is a propane chain ({parent_chain} carbons).")
    print(f"2. The amine group is located on carbon {amine_position}, so the parent name is propan-2-amine.")
    print(f"3. A phenyl group is attached as a substituent on carbon {phenyl_position}.")
    print(f"4. The complete IUPAC name is {phenyl_position}-phenylpropan-{amine_position}-amine.")
    
    final_name = "1-phenylpropan-2-amine"
    return final_name

# Execute the analysis
final_answer = analyze_compound_and_determine_iupac_name()

# Restore original stdout
sys.stdout = old_stdout
# Get the captured output
output_text = captured_output.getvalue()
# Print the detailed analysis
print(output_text)
# Print the final answer in the required format
print(f"<<<{final_answer}>>>")