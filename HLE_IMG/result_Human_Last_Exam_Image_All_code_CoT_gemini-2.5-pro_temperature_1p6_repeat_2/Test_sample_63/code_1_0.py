import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer to capture output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

print("Step-by-step analysis to identify the correct molecular structure from the NMR spectrum.")
print("-" * 70)

# Step 1: Analyze the observed 1H NMR Spectrum
print("Step 1: Analysis of the provided 1H NMR spectrum.")
print("The spectrum shows the following distinct signals:")
print(" - Signal at ~9.0 ppm: Broad singlet, integration ~1H. This is characteristic of an amide N-H proton.")
print(" - Signal at ~7.1 ppm: Multiplet, located in the aromatic region.")
print(" - Signal at ~3.4 ppm: Singlet, integration ~2H. Likely a CH2 group with no adjacent protons.")
print(" - Signal at ~2.8 ppm: Quartet, integration ~4H.")
print(" - Signal at ~2.3 ppm: Singlet, integration ~3H. Likely a CH3 group with no adjacent protons.")
print(" - Signal at ~1.1 ppm: Triplet, integration ~6H.")
print("\nKey Spectroscopic Pattern:")
print("The quartet at ~2.8 ppm (4H) and the triplet at ~1.1 ppm (6H) are a classic 'ethyl group' pattern.")
print("The n+1 rule confirms this: the CH2 group (4H) is split by the CH3 group (3 protons -> 3+1 = 4 lines, a quartet) and the CH3 group (6H) is split by the CH2 group (2 protons -> 2+1 = 3 lines, a triplet).")
print("This strongly indicates the presence of two equivalent ethyl groups attached to a nitrogen: a diethylamino group, -N(CH2CH3)2.")
print("-" * 70)


# Step 2: Evaluate the candidate structures
print("Step 2: Evaluating the candidate structures against the spectrum.")

print("\nAnalysis of Structures A-G and B-G:")
print(" - Both A-G and B-G contain a dimethylamino group, -N(CH3)2, not a diethylamino group.")
print(" - A -N(CH3)2 group would show up as a single peak (singlet) for 6 protons.")
print(" - The spectrum clearly shows the quartet-triplet pattern of a -N(C2H5)2 group.")
print("   => Conclusion: Structures A-G and B-G are incorrect.")

print("\nAnalysis of Structures C-L and D-L:")
print(" - Both C-L and D-L contain the required amide N-H, an aromatic ring, a -CO-CH2- group, and a diethylamino group, matching the general features.")
print(" - The key difference is the number of methyl groups on the aromatic ring.")
print("   - Structure C-L has one aromatic methyl group (-CH3).")
print("   - Structure D-L has two aromatic methyl groups (-CH3).")
print("\nLet's use the integration to distinguish between C-L and D-L:")
print(" - The signal at ~2.3 ppm is a singlet, characteristic of the aromatic methyl group(s).")
print(" - We can use the triplet at 1.1 ppm as a reference, which corresponds to 6 protons from the two ethyl -CH3 groups.")
print(" - The area of the singlet at 2.3 ppm is visually about half the area of the triplet at 1.1 ppm.")
print("   - This implies an integration ratio of 3 protons to 6 protons.")
print("   - An integration of 3H for the singlet at 2.3 ppm matches structure C-L, which has one aromatic methyl group (3H).")
print("   - Structure D-L has two aromatic methyl groups, which would result in a singlet with an integration of 6H. This is inconsistent with the spectrum.")
print("   => Conclusion: Structure C-L is the correct match, and D-L is incorrect.")
print("-" * 70)

# Step 3: Final Verification of Structure C-L
print("Step 3: Final verification of all signals for Structure C-L.")
print("Predicted signals for C-L vs. Observed signals:")
print(" - Amide NH:             Predicted: 1H, singlet.  Observed: ~9.0 ppm, singlet (1H).    => Match!")
print(" - Aromatic CH3:         Predicted: 3H, singlet.  Observed: ~2.3 ppm, singlet (3H).    => Match!")
print(" - -CO-CH2-N-:           Predicted: 2H, singlet.  Observed: ~3.4 ppm, singlet (2H).    => Match!")
print(" - -N-(CH2-CH3)2:        Predicted: 4H quartet + 6H triplet. Observed: ~2.8 ppm (q, 4H) + ~1.1 ppm (t, 6H). => Perfect Match!")
print(" - Aromatic H's:         Predicted: 4H, multiplet. Observed: ~7.1 ppm, multiplet. Integration is consistent. => Match!")
print("\nFinal Conclusion: The NMR spectrum unequivocally corresponds to structure C-L.")

# Get the captured output
final_output = output_buffer.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Close the buffer
output_buffer.close()

# Print the captured output
print(final_output)
<<<C>>>