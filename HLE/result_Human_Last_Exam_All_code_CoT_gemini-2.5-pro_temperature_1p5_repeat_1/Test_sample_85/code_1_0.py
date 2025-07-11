import sys

# This script is a "literate program" that explains the reasoning for the linguistic reconstruction.
# It doesn't compute the answer algorithmically but rather demonstrates the logic step-by-step.

def solve_reconstruction():
    """
    Reconstructs the Sanskrit ancestor for two invented Indo-Aryan words.
    """
    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"

    # --- Reconstruction Components ---
    reconstructed_initial = "k"
    reconstructed_vowel = "a"
    reconstructed_cluster = "nśy"
    
    # --- Final Proto-form Assembly ---
    sanskrit_stem = reconstructed_initial + reconstructed_vowel + reconstructed_cluster
    sanskrit_ipa = "kanɕja"

    print("Reconstruction of the Sanskrit ancestor for:")
    print(f"  - Hindi: {hindi_word}")
    print(f"  - Punjabi: {punjabi_word}\n")

    print("--- Step 1: Initial Consonant ---")
    print("The initial consonant in both words is 'k'. This sound is stable from Sanskrit.")
    print(f"Reconstructed initial: *{reconstructed_initial}\n")

    print("--- Step 2: Vowel and Nasalization ---")
    print("Hindi 'ãː' (long) and Punjabi 'ə̃' (short) both show nasalization ('̃'), pointing to a nasal in the original consonant cluster.")
    print("The long 'aː' in Hindi results from compensatory lengthening, which occurs when a consonant cluster simplifies. This implies the original Sanskrit vowel was short.")
    print(f"Reconstructed vowel: *{reconstructed_vowel}\n")

    print("--- Step 3: The Consonant Cluster ---")
    print("The central puzzle is the cluster that yields '-s' in Hindi and '-d͡ʒʱ' in Punjabi.")
    print("A plausible Sanskrit cluster that explains this divergence is *nśy.")
    print(" - Hindi path: *nśy > Pkt. *ṃss > Hindi *s (with lengthening of the preceding vowel).")
    print(" - Punjabi path: *nśy undergoes palatal assimilation, becoming voiced and aspirated (e.g. > *ñjh > *d͡ʒʱ).\n")

    print("--- Step 4: Final Assembled Form (The Equation) ---")
    print("Combining the reconstructed parts gives the final proto-word stem:")
    print(f"  Initial ('{reconstructed_initial}') + Vowel ('{reconstructed_vowel}') + Cluster ('{reconstructed_cluster}') => Proto-stem (*{sanskrit_stem})\n")
    
    print("--- Final Answer ---")
    print("The most probable Sanskrit ancestor, represented as a stem in IPA, is:")
    print(sanskrit_ipa)


# Execute the function to print the solution
solve_reconstruction()

# The final answer in the requested format will follow the code block.
# This ensures the output is a single string in IPA as requested.
final_answer_ipa = "kanɕja"
sys.stdout.write(f"\n<<<kancja>>>\n")