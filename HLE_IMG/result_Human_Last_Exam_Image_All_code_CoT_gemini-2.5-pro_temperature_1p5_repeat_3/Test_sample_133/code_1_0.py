import sys
import io

# Helper function to print the output
def print_analysis():
    """
    Performs and prints the step-by-step analysis of the mass spectrum.
    """
    print("Step 1: Analysis of the Molecular Ion Peak")
    print("------------------------------------------")
    m_z_262 = 29.0
    m_z_264 = 12.9
    m_z_266 = 3.22
    
    print("The peaks at the highest mass-to-charge ratio (m/z) are observed at 262, 264, and 266.")
    print(f"Their relative intensities are approximately I(262)={m_z_262}%, I(264)={m_z_264}%, I(266)={m_z_266}%.")
    print("The presence of M, M+2, and M+4 peaks strongly suggests the presence of two chlorine (Cl) atoms, due to the natural isotopic abundance of ³⁵Cl (~75.8%) and ³⁷Cl (~24.2%).")
    print("The theoretical intensity ratio for a fragment with two chlorines (Cl₂) is approximately 9:6:1.")
    print("While the observed ratio is not a perfect match, this pattern is the most plausible explanation for the cluster.")
    print("Therefore, we assume the molecular ion peak M⁺ (containing two ³⁵Cl isotopes) is at m/z = 262.")
    
    mw = 262
    cl_mass = 2 * 35
    remainder_mass = mw - cl_mass
    print(f"\nThe molecular weight of the compound (with the most abundant isotopes) is {mw} u.")
    print(f"The mass contribution from two ³⁵Cl atoms is 2 * 35 = {cl_mass} u.")
    print(f"The remaining mass for the CₓHᵧ part of the molecule is {mw} - {cl_mass} = {remainder_mass} u.")
    print("\nA possible molecular formula that fits this mass is C₁₅H₁₂Cl₂ (15*12 + 12*1 + 2*35 = 180 + 12 + 70 = 262).")
    
    print("\nStep 2: Identification of the Base Peak")
    print("--------------------------------------")
    base_peak_mz = 225
    print(f"The base peak (the most intense peak, I=100%) is at m/z = {base_peak_mz}.")
    
    print("\nStep 3: Analysis of Fragmentation Pathways")
    print("------------------------------------------")
    fragment_227_mz = 227
    loss_from_M_to_227 = mw - fragment_227_mz
    print(f"A major fragment is seen at m/z = {fragment_227_mz}. The loss from the molecular ion is {mw} - {fragment_227_mz} = {loss_from_M_to_227} u.")
    print("This corresponds to the loss of a chlorine radical (³⁵Cl), a very common fragmentation for chlorinated compounds. This creates the [M-Cl]⁺ ion.")
    print("The corresponding isotope peak for [M-Cl]⁺ containing one ³⁷Cl is found at m/z = 229 (Intensity = 22.2%). The parent peak at m/z=227 has intensity 66.7%. The ratio 22.2/66.7 ≈ 0.33 matches the expected ³⁷Cl/³⁵Cl ratio of ~1:3 perfectly.")
    
    loss_from_227_to_225 = fragment_227_mz - base_peak_mz
    print(f"\nThe base peak at m/z = {base_peak_mz} is formed from the fragment at m/z = {fragment_227_mz} by a loss of {fragment_227_mz} - {base_peak_mz} = {loss_from_227_to_225} u.")
    print("This corresponds to the loss of a hydrogen molecule (H₂), which is plausible as it can lead to a more stable, rearranged ion.")

    fragment_190_mz = 190
    loss_from_225_to_190 = base_peak_mz - fragment_190_mz
    print(f"\nAnother significant fragment at m/z = {fragment_190_mz} can be explained by a loss from the base peak. The mass difference is {base_peak_mz} - {fragment_190_mz} = {loss_from_225_to_190} u.")
    print("This corresponds to the loss of the second chlorine atom (³⁵Cl) from the base peak fragment.")
    
    print("\nSummary of the fragmentation pathway:")
    print(f"C₁₅H₁₂Cl₂⁺ (m/z = {mw})  -> [loss of Cl] -> C₁₅H₁₂Cl⁺ (m/z = {fragment_227_mz}) -> [loss of H₂] -> C₁₅H₁₀Cl⁺ (m/z = {base_peak_mz}) -> [loss of Cl] -> C₁₅H₁₀⁺ (m/z = {fragment_190_mz})")

    print("\nStep 4: Proposed Structure and IUPAC Name")
    print("------------------------------------------")
    print("The formula C₁₅H₁₂Cl₂ and the fragmentation pattern are consistent with a bis(chlorophenyl)cyclopropane structure.")
    print("Given the symmetrical fragmentation, a likely isomer is 1,1-bis(4-chlorophenyl)cyclopropane.")
    
    final_name = "1,1-bis(4-chlorophenyl)cyclopropane"
    print(f"\nThe compound corresponding to the mass spectrum is: {final_name}")

# Execute the analysis
print_analysis()

# Final Answer Block
final_answer = "1,1-bis(4-chlorophenyl)cyclopropane"
# Suppressing direct output to screen for the final answer block as per user request to use 'print'
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()
print(f"<<<{final_answer}>>>")
sys.stdout = old_stdout
# The capture is just for adhering to the "no direct output" instruction, the actual value is in the variable.
# For the purpose of this simulation, we'll print what would be the final line.
print("\n<<<1,1-bis(4-chlorophenyl)cyclopropane>>>")