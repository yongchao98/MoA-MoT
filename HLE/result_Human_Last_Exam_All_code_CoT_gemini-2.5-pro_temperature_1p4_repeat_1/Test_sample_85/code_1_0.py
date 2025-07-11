import sys
# Set default encoding to UTF-8 to ensure IPA characters print correctly.
sys.stdout.reconfigure(encoding='utf-8')

def reconstruct_sanskrit_ancestor():
    """
    This function explains and prints the reconstruction of a Sanskrit word
    from its invented Hindi and Punjabi descendants.
    """
    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"
    proto_sanskrit_word = "kaṁkṣa"
    proto_ipa = "kaŋkʂa"

    print(f"Reconstruction of Sanskrit Ancestor for Hindi '{hindi_word}' and Punjabi '{punjabi_word}'")
    print("-" * 75)
    print(f"Proposed Sanskrit Proto-form: {proto_sanskrit_word}")
    print(f"This form means 'desire' or 'longing' in Sanskrit.")

    print("\nDerivation Path to Hindi 'kãːs':")
    print(f"1. Sanskrit: {proto_sanskrit_word} (/kaŋkʂa/)")
    print("   The Sanskrit cluster 'kṣ' has several different outcomes in its daughter languages.")
    print("2. Prakrit (Western Dialect Branch): > kaṁska")
    print("   In this branch, 'kṣ' becomes 'sk'.")
    print("3. Late Prakrit / Old Hindi: > kaṁsa")
    print("   The consonant cluster 'sk' simplifies to 's'.")
    print("4. Modern Hindi: > kãːs")
    print("   The preceding vowel 'a' lengthens to 'aː' to compensate for the cluster simplification. The nasal 'ṁ' is realized as nasalization on the vowel. The final short 'a' is dropped.")

    print("\nDerivation Path to Punjabi 'kə̃d͡ʒʱ':")
    print(f"1. Sanskrit: {proto_sanskrit_word} (/kaŋkʂa/)")
    print("   The same Sanskrit cluster 'kṣ' follows a different path in another dialect.")
    print("2. Prakrit (Minority Dialectal Form): > kaṁjha")
    print("   A documented, though less common, change is 'kṣ' becoming 'jh'.")
    print("3. Old Punjabi: > kəñjh")
    print("   The vowel 'a' shortens to 'ə', a common change in Punjabi before consonant clusters.")
    print("4. Modern Punjabi: > kə̃d͡ʒʱ")
    print("   The palatal aspirate 'jh' strengthens to a voiced aspirated affricate 'd͡ʒʱ'. The pre-consonantal nasal becomes nasalization on the vowel.")

    print("-" * 75)
    print("Conclusion: The single proto-form 'kaṁkṣa' can account for all phonological features in both descendants via divergent but attested sound changes.")
    print("\nFinal Reconstructed Form in IPA:")
    print(proto_ipa)


reconstruct_sanskrit_ancestor()
