import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

print("Analysis of Spectroscopic Data to Determine IUPAC Name")
print("="*60)

# 1. Mass Spectrometry (MS) Analysis
print("\nStep 1: Mass Spectrum Analysis")
print(" - Molecular ion (M+) peak is observed at m/z = 135.")
print(" - An odd molecular weight suggests the presence of an odd number of Nitrogen atoms (Nitrogen Rule).")
mw = 135
print(f" - Assuming 1 Nitrogen atom.")

# 2. Molecular Formula and Degree of Unsaturation
print("\nStep 2: Molecular Formula and Degree of Unsaturation")
# From full analysis: 9 carbons, 13 hydrogens, 1 nitrogen
C = 9
H = 13
N = 1
calculated_mw = C*12 + H*1 + N*14
print(f" - Based on all spectra, the molecular formula is deduced to be C{C}H{H}N.")
print(f" - Calculated MW for C{C}H{H}N = ({C}*12) + ({H}*1) + ({N}*14) = {calculated_mw}, which matches the m/z of {mw}.")
unsaturation = C - H/2 - N/2 + 1
print(f" - Degree of Unsaturation = C - H/2 - N/2 + 1 = {C} - {H}/2 - {N}/2 + 1 = {unsaturation}.")
print(" - A value of 4.5 suggests a calculation error, let's recheck the formula. The correct formula is (2C+2+N-H)/2.")
unsaturation_correct = (2*C + 2 + N - H) / 2
print(f" - Correct Degree of Unsaturation = (2*C + 2 + N - H) / 2 = (2*{C} + 2 + {N} - {H}) / 2 = {int(unsaturation_correct)}.")
print(" - A degree of unsaturation of 4 suggests a benzene ring (4 degrees: 3 double bonds + 1 ring).")

# 3. IR and NMR Analysis
print("\nStep 3: IR and NMR Analysis")
print(" - IR spectrum: Shows N-H stretches (~3300-3400 cm-1, primary amine), aromatic C-H (>3000 cm-1), and aliphatic C-H (<3000 cm-1).")
print(" - 1H NMR: Shows a monosubstituted phenyl group (~7.3 ppm, 5H), a CH3 doublet (~1.2 ppm, 3H), and two aliphatic multiplets (1H and 2H).")
print(" - 13C NMR & DEPT-135: Indicate 7 distinct carbon environments.")
c_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
print(f"   - Signals: {c_shifts}")
print("   - Aromatic region (126-146 ppm): 4 signals representing 6 carbons of a monosubstituted benzene ring.")
print("   - Aliphatic region (19-50 ppm): 3 signals.")
print("   - DEPT-135 shows one negative signal (CH2), and five positive signals (CH, CH3).")
print("   - HSQC correlation confirms: CH3 at 19.2 ppm, CH2 at 43.5 ppm, and CH at 49.6 ppm.")

# 4. Structure Elucidation
print("\nStep 4: Structure Elucidation")
print(" - The data indicates a phenyl group and a three-carbon chain with a primary amine group.")
print(" - The fragment -CH-CH3 is confirmed by the 3H doublet in the 1H NMR.")
print(" - The placement of the -NH2 and C6H5 groups leads to two main isomers:")
print("   A) 1-phenylpropan-2-amine: C6H5-CH2-CH(NH2)-CH3")
print("   B) 2-phenylpropan-1-amine: C6H5-CH(CH3)-CH2-NH2")
print(" - In structure A, the CH carbon is attached to the electronegative NH2 group, and the CH2 is benzylic.")
print(" - This matches the 13C data, where the CH carbon (49.6 ppm) is more deshielded (further downfield) than the CH2 carbon (43.5 ppm).")
print(" - Therefore, the correct structure is 1-phenylpropan-2-amine.")

# 5. IUPAC Naming
print("\nStep 5: IUPAC Naming")
print(" - Parent chain: propane (3 carbons).")
print(" - Substituents: an amine group on carbon 2 and a phenyl group on carbon 1.")
print(" - Final Name: 1-phenylpropan-2-amine")
print("="*60)

# Restore stdout
sys.stdout = old_stdout
# Print the captured output to the actual console
output_string = captured_output.getvalue()
print(output_string)

# Final answer format
print("<<<1-phenylpropan-2-amine>>>")