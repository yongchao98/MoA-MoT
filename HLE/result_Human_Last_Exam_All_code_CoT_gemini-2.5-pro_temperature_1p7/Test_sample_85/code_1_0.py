def reconstruct_sanskrit_from_indo_aryan(hindi_word, punjabi_word):
    """
    Reconstructs and explains the Sanskrit proto-word for a given set of
    cognates from Hindi and Punjabi.

    This function applies principles of historical linguistics to deduce the
    most probable ancestor form.
    """
    print(f"Task: Reconstruct the Sanskrit ancestor of Hindi '{hindi_word}' and Punjabi '{punjabi_word}'.\n")

    # The hypothesized Sanskrit ancestor in IAST (International Alphabet of Sanskrit Transliteration)
    # and IPA (International Phonetic Alphabet). The hypothesis is based on analyzing the
    # divergent outcomes of a specific consonant cluster.
    proto_form_iast = "kaṁsya"
    proto_form_ipa = "kãsja"

    print("--- Step 1: Analysis of Sound Correspondences ---")
    print("  - Initial 'k' is stable, pointing to Sanskrit 'k'.")
    print("  - Vowel nasality ('ãː', 'ə̃') points to a Sanskrit nasal (e.g., 'ṁ', 'n', 'ñ').")
    print("  - The correspondence of Hindi 's' and Punjabi 'd͡ʒʱ' is the key.")
    print("  - This suggests a complex cluster, likely '-ṁsya-', which has different reflexes.\n")

    print("--- Step 2: Derivation Path to Hindi ---")
    print(f"1. Proposed Sanskrit Form: {proto_form_iast} (IPA: {proto_form_ipa})")
    print("2. Development into Middle Indo-Aryan (Prakrit): *kaṁssa")
    print("   (Reason: The cluster 'sy' regularly becomes a geminate 'ss' in Prakrit.)")
    print(f"3. Development into Hindi: {hindi_word}")
    print("   (Reason: The geminate 'ss' simplifies to 's', causing compensatory lengthening of the preceding vowel 'a' to 'aː'. The nasal 'ṁ' is retained as vowel nasality.)\n")

    print("--- Step 3: Derivation Path to Punjabi ---")
    print(f"1. Proposed Sanskrit Form: {proto_form_iast} (IPA: {proto_form_ipa})")
    print("2. Development into Proto-Punjabi: *kañjha")
    print("   (Reason: A specific regional sound change, 'sya' > 'jha', is attested in western dialects. The nasal 'ṁ' assimilates to 'ñ' before the palatal 'jh'.)")
    print(f"3. Development into Punjabi: {punjabi_word}")
    print("   (Reason: The short vowel 'a' centralizes to 'ə'. The nasal 'ñ' becomes vowel nasality. The voiced aspirate 'jh' is preserved as 'd͡ʒʱ'.)\n")
    
    print("--- Step 4: Conclusion ---")
    print("Based on these derivations, the most probable Sanskrit ancestor is identified.")
    print("\nFinal Reconstructed Sanskrit Form in IPA:")
    print(proto_form_ipa)


# The input words from the problem description
hindi_descendant = "kãːs"
punjabi_descendant = "kə̃d͡ʒʱ"

reconstruct_sanskrit_from_indo_aryan(hindi_descendant, punjabi_descendant)