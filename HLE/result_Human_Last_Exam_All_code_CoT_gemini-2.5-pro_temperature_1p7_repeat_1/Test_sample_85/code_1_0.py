import pandas as pd

def reconstruct_sanskrit():
    """
    This function explains the step-by-step reconstruction of a Sanskrit word
    from its invented descendants in Hindi and Punjabi.
    """
    
    # Define the data in a structured way
    data = {
        'Language': ['Hindi', 'Punjabi', 'Proto-Sanskrit (Reconstruction)'],
        'IPA': ['kãːs', 'kə̃d͡ʒʱ', ''],
        'Notes': [
            'Descendant 1 (Central Indo-Aryan)', 
            'Descendant 2 (Northwestern Indo-Aryan)',
            'The form we want to reconstruct'
        ]
    }
    df = pd.DataFrame(data)
    
    print("--- Word Analysis ---")
    print(df[['Language', 'IPA', 'Notes']])
    print("\n--- Reconstruction Steps ---")
    
    print("1. Initial Consonant 'k':")
    print("   - Hindi has 'k'.")
    print("   - Punjabi has 'k'.")
    print("   - This is a direct correspondence. The proto-form must begin with *k.")
    
    print("\n2. Vowel and Nasalization ('ãː' vs 'ə̃'):")
    print("   - Both vowels are nasalized ('~'). This indicates the loss of a nasal consonant in the original word, with its nasal quality being transferred to the vowel.")
    print("   - The vowel length difference (Hindi 'aː' vs. Punjabi 'ə') is a regular correspondence. Sanskrit long 'ā' often remains long in Hindi but is shortened and centralized to 'ə' in Punjabi.")
    print("   - Therefore, we reconstruct a long vowel *ā followed by a nasal consonant.")

    print("\n3. Medial Consonant Cluster ('s' vs 'd͡ʒʱ'):")
    print("   - This is the most complex part. Hindi's 's' suggests simplification from a sibilant cluster. Punjabi's 'd͡ʒʱ' is a voiced, aspirated palatal affricate.")
    print("   - The palatal nature of Punjabi 'd͡ʒʱ' points to a palatal cluster in Sanskrit. The best candidate is 'ñś' ([ɲɕ]).")
    print("   - This cluster neatly explains the divergence:")
    print("     - Path to Hindi: The palatal cluster 'ñś' simplifies to 's' in the central languages. The palatal nasal 'ñ' is the source of the vowel nasalization. E.g., *kāñśa > *kā̃sa > kãːs.")
    print("     - Path to Punjabi: The cluster 'ñś' undergoes a special development in the northwestern branch, becoming 'd͡ʒʱ'. The palatal articulation is preserved, but the cluster is strengthened (fortition). The nasality of 'ñ' transfers to the preceding vowel. E.g., *kāñśa > kə̃d͡ʒʱ.")

    print("\n4. Final Reconstruction:")
    print("   - Combining the reconstructed parts: *k + *ā + *ñś + a common final vowel '-a'.")
    print("   - Proto-form (IAST): kāñśa")
    print("   - Proto-form (IPA): kaːɲɕa")

    # The reconstructed word in IPA
    final_reconstruction = "kaːɲɕa"
    print(f"\nFinal reconstructed Sanskrit ancestor in IPA is: {final_reconstruction}")

reconstruct_sanskrit()
<<<kaːɲɕa>>>