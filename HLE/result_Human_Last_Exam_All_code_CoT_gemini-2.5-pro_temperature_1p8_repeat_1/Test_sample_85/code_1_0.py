def reconstruct_sanskrit():
    """
    This function explains the reconstruction of a Sanskrit word from its
    hypothetical Hindi and Punjabi descendants and prints the result.
    """
    hindi_word = "kãːs"
    punjabi_word = "kənd͡ʒʱ"
    
    # Reconstructed components based on linguistic rules
    initial = "k"
    vowel1 = "a"
    nasal = "m"
    sibilant = "ɕ" # IPA for Sanskrit 'ś'
    vowel2 = "a"
    
    reconstructed_word = f"/{initial}{vowel1}{nasal}{sibilant}{vowel2}/"

    print("Reconstruction of the Sanskrit Ancestor")
    print("-" * 40)
    print(f"Descendants:")
    print(f"  - Hindi: {hindi_word}")
    print(f"  - Punjabi: {punjabi_word}")
    print("\nReasoning:")
    print("1. Initial 'k' is preserved in both languages.")
    print("2. The correspondence of Hindi '-ā̃s' and Punjabi '-ə̃jh' (IPA -ənd͡ʒʱ) points to the Sanskrit cluster '-aṁś-'.")
    print("   - The nasalization ('̃') comes from the Sanskrit anusvāra 'ṁ'.")
    print("   - The sibilant 's' in Hindi and the affricate 'jh' (d͡ʒʱ) in Punjabi both derive from Sanskrit 'ś' in this context.")
    print("   - The long 'ā' in Hindi is compensatory lengthening from the cluster simplification.")
    print("3. The final short '-a' of the Sanskrit form is dropped in both modern languages, a common feature.")
    print("\nReconstructed Proto-Word (in IPA):")
    print(f"The reconstructed components are: /{initial}/ + /{vowel1}/ + /{nasal}/ + /{sibilant}/ + /{vowel2}/")
    print(f"\nFinal Answer: {reconstructed_word}")

reconstruct_sanskrit()
<<< /kamɕa/ >>>