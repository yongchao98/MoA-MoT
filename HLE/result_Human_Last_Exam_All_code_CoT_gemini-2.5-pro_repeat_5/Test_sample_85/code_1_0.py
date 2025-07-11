import unicodedata

def normalize_ipa(ipa_string):
    """Normalize IPA string to handle different representations."""
    # Using NFD to decompose characters like '̃' into base + combining mark
    return unicodedata.normalize('NFD', ipa_string)

def reconstruct_sanskrit_ancestor():
    """
    This function reconstructs the Sanskrit ancestor for Hindi kãːs and Punjabi kə̃d͡ʒʱ.
    The process involves analyzing sound correspondences and applying historical linguistic principles.
    """
    
    # Descendant words in IPA
    hindi_word_ipa = "kãːs"
    punjabi_word_ipa = "kə̃d͡ʒʱ"

    # Step 1: Initial consonant correspondence
    # Hindi k -> Punjabi k. Proto-form likely starts with 'k'.
    reconstruction_step1 = "k"

    # Step 2: Vowel and nasalization correspondence
    # Hindi ãː -> Punjabi ə̃. Suggests a Sanskrit long 'ā' with a preceding or following nasal.
    # The length is preserved in Hindi, shortened in Punjabi. Nasalization comes from 'ṁ'.
    reconstruction_step2 = "kaːm" # Using 'm' as a common representation for 'ṁ' before the cluster.

    # Step 3: Medial cluster correspondence
    # Hindi s <-> Punjabi d͡ʒʱ.
    # The most plausible source cluster is 'ṁś' (nasal + palatal sibilant).
    # Path to Hindi: ṁś > ṁs > ̃s (with compensatory lengthening, though 'ā' is already long).
    # Path to Punjabi: ṁś > ̃d͡ʒʱ (irregular change involving voicing, affrication, and aspiration).
    # This leads to a proto-form like 'kāṁśa'.
    reconstruction_step3 = "kaːmɕa" # Using IPA for śa.

    # The reconstructed word in Devanagari is 'कांश्य' or more likely 'कांष' for this problem.
    # A better representation is kāṁśa.
    # The IPA for kāṁśa is /kaːɲɕə/. The nasal 'ṁ' becomes palatal 'ɲ' before the palatal 'ɕ'.
    # The final 'a' in Sanskrit is pronounced as a schwa 'ə'.
    final_reconstruction_ipa = "kaːɲɕə"

    print(f"Hindi descendant: {hindi_word_ipa}")
    print(f"Punjabi descendant: {punjabi_word_ipa}")
    print("-" * 20)
    print("Reconstruction Steps:")
    print("1. Initial 'k' is stable in both languages -> Proto-form starts with 'k'.")
    print("2. Vowel correspondence (Hindi ãː vs. Punjabi ə̃) points to a Sanskrit long 'ā' with a nasal -> Proto-form has 'kāṁ...'.")
    print("3. Medial consonant correspondence (Hindi 's' vs. Punjabi 'd͡ʒʱ') is best explained by a proto-cluster 'ṁś'.")
    print("   - For Hindi: Skt. ṁś > Pkt. ṁs > Hindi ̃s. (e.g., 'vaṁśa' > 'bā̃s'). This is a regular sound change.")
    print("   - For Punjabi: Skt. ṁś > Pbi. ̃d͡ʒʱ. This is an irregular change, but more plausible than alternatives.")
    print("4. The full reconstructed Sanskrit word is 'kāṁśa'.")
    print("-" * 20)
    print("Final reconstructed ancestor in IPA:")
    # The final output needs to be just the IPA string per the instructions, but enclosed in the special format.
    # To satisfy the "output each number in the final equation" rule in a creative way,
    # let's print the characters of the final IPA string.
    print(f"The final reconstruction is composed of the following IPA characters:")
    for char in final_reconstruction_ipa:
        print(char)
    
    # Final answer in the required format
    print(f"\nFinal Answer String: {final_reconstruction_ipa}")

reconstruct_sanskrit_ancestor()
print("<<<kaːɲɕə>>>")