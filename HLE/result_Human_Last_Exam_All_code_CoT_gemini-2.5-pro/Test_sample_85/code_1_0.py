import unicodedata

def solve_reconstruction():
    """
    This function explains the process of reconstructing the Sanskrit ancestor
    and prints the final result in IPA.
    """
    
    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"
    
    # Step 1: Analyze the initial consonant
    initial_consonant = "k"
    
    # Step 2: Analyze the vowel and nasality
    # The correspondence of Hindi ā̃ vs Punjabi ə̃ points to an original short 'a' before a consonant cluster that contained a nasal.
    vowel = "a"
    
    # Step 3: Analyze the final consonant cluster
    # The correspondence of Hindi 's' vs Punjabi 'd͡ʒʱ' is the main puzzle.
    # It points to a complex Sanskrit cluster with divergent reflexes. The most plausible candidate is 'kṣ' (/kʂ/).
    # The nasal that precedes a 'k' in Sanskrit is the velar nasal 'ṅ' (/ŋ/).
    cluster = "ṅkṣ" # Spelled 'nksh' in IAST, /ŋkʂ/ in IPA
    
    # Step 4: Combine the elements to form the Sanskrit ancestor
    # The reconstructed form is kaṅkṣa, which includes a final short 'a' common in Sanskrit nouns.
    reconstructed_word_devanagari = "कङ्क्ष"
    reconstructed_word_iast = "kaṅkṣa"
    
    # Step 5: Convert the final reconstruction to IPA
    # k -> k
    # a -> a
    # ṅ -> ŋ
    # kṣ -> kʂ
    # a -> a
    # Note: For IPA representation, the final short 'a' in Sanskrit is often transcribed as /ɐ/. We will use 'a' for simplicity as it's phonemic.
    final_ipa = "kaŋkʂa"

    print("Reconstruction Steps:")
    print(f"1. Initial consonant: The 'k' in '{hindi_word}' and '{punjabi_word}' points to an ancestor with 'k'.")
    print(f"2. Vowel and Nasality: The 'ãː'/'ə̃' correspondence indicates an original short 'a' before a nasal + consonant cluster.")
    print(f"3. Final Consonant Cluster: The 's'/'d͡ʒʱ' correspondence is best explained by the complex Sanskrit cluster 'kṣ', preceded by the nasal 'ṅ'.")
    print(f"4. Reconstructed Sanskrit form: {reconstructed_word_iast}")
    print("\nFinal Answer in IPA:")
    print(final_ipa)
    
    # The required final output format
    print(f"\n<<<{final_ipa}>>>")

solve_reconstruction()