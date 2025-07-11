import io
import sys

# Mock sys.stdout to capture printed output
stdout_capture = io.StringIO()
sys.stdout = stdout_capture

# This script reconstructs the Sanskrit ancestor and prints the reasoning.

# --- Step 1: Analyze the Descendant Words ---
hindi_word = "kãːs"
punjabi_word = "kə̃d͡ʒʱ"

print("Reconstruction Plan:")
print("1. Compare the initial consonants.")
print("2. Compare the vowels and nasalization.")
print("3. Compare the medial consonants/clusters to find the most likely source.")
print("4. Combine the reconstructed sounds into a single Sanskrit word in IPA.\n")

# --- Step 2: Step-by-step Reconstruction ---
print("Executing the Reconstruction:")

# Initial consonant
initial_k = "k"
print(f"-> The initial consonant in both Hindi ({hindi_word}) and Punjabi ({punjabi_word}) is '{initial_k}'.")
print(f"   Therefore, the reconstructed initial is *{initial_k}.")

# Vowel and Nasalization
sanskrit_vowel_nasal = "aṃ"
print(f"-> The Hindi 'ãː' and Punjabi 'ə̃' both point to an original short 'a' plus a nasal, which we reconstruct as *{sanskrit_vowel_nasal}.")

# Medial Consonant Cluster
sanskrit_consonant = "ś"
print(f"-> The Hindi '-s' and Punjabi '-d͡ʒʱ' correspondence is key.")
print(f"   The sound change *ṃś > -̃s is regular in Hindi.")
print(f"   The sound change *ṃś > -̃d͡ʒ(ʱ) is known in Punjabi (e.g., Sanskrit vaṃśa > Punjabi və̃d͡ʒ).")
print(f"   This points to a Sanskrit palatal sibilant *{sanskrit_consonant} following the nasal.")

# --- Step 3: Final Reconstruction ---
reconstructed_stem = "kaṃś"
reconstructed_word = "kaṃśa"
# IPA Transcription: k + ɐ̃ (nasalized a) + ɕ (palatal sibilant 'ś') + ɐ (final vowel 'a')
reconstructed_ipa = "kɐ̃ɕɐ"

print(f"\nConclusion:")
print(f"Combining the parts ({initial_k} + {sanskrit_vowel_nasal} + {sanskrit_consonant}), we get the stem *{reconstructed_stem}.")
print(f"The most probable Sanskrit ancestor word is '{reconstructed_word}'.")
print(f"In IPA, the final reconstructed form is: {reconstructed_ipa}")

# Restore original stdout and print captured output to the user
sys.stdout = sys.__stdout__
output = stdout_capture.getvalue()
print(output)

# Final answer in the required format
# This print statement is for the platform to capture the final answer.
# It will not be seen by the user in the final output block.
final_answer_string_for_capture = f"<<<{reconstructed_ipa}>>>"
# The line above is a special instruction for the system.
# The final formatted answer will be extracted from it.