def reconstruct_sanskrit_ancestor():
    """
    Prints the step-by-step reconstruction of a Sanskrit word from its
    hypothetical Indo-Aryan descendants.
    """

    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"
    
    print("Reconstruction of the Sanskrit ancestor for Hindi '{}' and Punjabi '{}'".format(hindi_word, punjabi_word))
    print("-" * 70)
    
    print("Step 1: Analyze the sound correspondences between the descendants.")
    print("  - Initial Consonant: Hindi 'k' corresponds to Punjabi 'k'.")
    print("  - Vowel & Nasalization: Hindi 'ãː' (long, nasal) corresponds to Punjabi 'ə̃' (short, nasal).")
    print("  - Medial Consonant: Hindi 's' corresponds to Punjabi 'd͡ʒʱ'.\n")
    
    print("Step 2: Reconstruct the Sanskrit sounds based on these correspondences.")
    print("  - Initial 'k': The correspondence is direct. The Sanskrit ancestor began with 'k'.")
    print("  - Vowel/Nasal: The long 'ãː' in Hindi suggests compensatory lengthening from an original short '*a' + Nasal cluster. The Punjabi 'ə̃' supports an original short '*a' + Nasal. So, we reconstruct '*aN...' (where N is a nasal).\n")

    print("Step 3: Reconstruct the medial consonant cluster. This is the most complex part.")
    print("  - The correspondence 's' ~ 'd͡ʒʱ' is unusual but points to a single Sanskrit cluster with two different historical developments (reflexes).")
    print("  - The Sanskrit cluster 'kṣ' (/kʂ/) is a known source for such divergent reflexes.")
    print("    - Path 1 (to Hindi): 'kṣ' is known to simplify to 's' in some contexts (e.g., Skt. kakṣa > 'kā̃sī' - a type of grass/scythe).")
    print("    - Path 2 (to Punjabi): 'kṣ' is also known to develop into 'jh' (e.g., Skt. kṣarati > Pkt. jharai). This 'jh' regularly becomes 'd͡ʒʱ' in Punjabi.")
    print("  - Therefore, the Sanskrit cluster was most likely '*kṣ'.\n")

    print("Step 4: Assemble the full Sanskrit word.")
    print("  - Combining the parts: k + a + Nasal + kṣ + a (standard stem ending).")
    print("  - In Sanskrit, a nasal before a velar stop like 'k' is written as the velar nasal 'ṅ'.")
    print("  - The reconstructed word is therefore: kaṅkṣa.\n")

    print("Step 5: Provide the final answer in IPA.")
    print("  - The Sanskrit word 'kaṅkṣa' is transcribed phonetically into IPA.")
    
    sanskrit_reconstruction = "kaṅkṣa"
    ipa_transcription = "/kaŋkʂa/"
    
    print("-" * 70)
    print("Final Reconstructed Form (Sanskrit): {}".format(sanskrit_reconstruction))
    print("Final Reconstructed Form (IPA): {}".format(ipa_transcription))
    print("-" * 70)

if __name__ == '__main__':
    reconstruct_sanskrit_ancestor()
