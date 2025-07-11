import sys

def reconstruct_sanskrit():
    """
    This function outlines the step-by-step reconstruction of a Sanskrit word
    from its invented Hindi and Punjabi descendants and prints the final result in IPA.
    """

    # --- Introduction ---
    print("Reconstructing the Sanskrit ancestor for Hindi 'kãːs' and Punjabi 'kə̃d͡ʒʱ'.")
    print("The reconstruction is based on comparing the phonological components of the two words.\n")

    # --- Step-by-step breakdown ---
    print("--- Reconstruction Steps ---")

    # Step 1: Initial Consonant
    print("1. Reconstructing the Initial Sound:")
    part_1 = "*k"
    print(f"   - Hindi 'k' and Punjabi 'k' both point to a stable initial sound.")
    print(f"   - Reconstructed Initial: {part_1}\n")

    # Step 2: Vowel and Nasalization
    print("2. Reconstructing the Vowel:")
    part_2 = "*ãː"
    print(f"   - Hindi 'ãː' (long, nasal) and Punjabi 'ə̃' (short, nasal) suggest a Sanskrit long vowel 'ā' followed by an anusvāra 'ṃ'.")
    print(f"   - The difference in vowel length is explained by the simplification of a following consonant cluster in Hindi, which did not happen in Punjabi.")
    print(f"   - Reconstructed Vowel Element: *āṃ (pronounced roughly as {part_2})\n")

    # Step 3: Medial Consonant Cluster
    print("3. Reconstructing the Medial Cluster:")
    part_3 = "*sja"
    print(f"   - Hindi 's' and Punjabi 'd͡ʒʱ' are reflexes of the Sanskrit cluster '-sya'.")
    print(f"   - This specific sound change is confirmed by the word for 'bronze': Skt. 'kāṃsya' -> Hi. 'kã̄sā', Pj. 'kañjh'.")
    print(f"   - Reconstructed Cluster: {part_3}\n")
    
    # Step 4: Final Assembled Word in IPA
    print("--- Final Reconstructed Form in IPA ---")
    final_reconstruction = "kãːsja"
    print("Combining the reconstructed parts (*k + *āṃ + *sya) gives us the Sanskrit word 'kāṃsya'.")
    print("In the International Phonetic Alphabet (IPA), 'kāṃsya' is transcribed as:")
    
    # Final equation/breakdown of the resulting IPA string
    # We "output each number [part] in the final equation"
    print(f"k + ãː + s + j + a => {final_reconstruction}")
    

reconstruct_sanskrit()

# Suppress the default print behavior after running the script if not in an interactive shell
# This is a good practice for clean output in some environments.
# In this case, we just need to provide the final answer tag.
sys.stdout.flush()

print("\n<<<kãːsja>>>")