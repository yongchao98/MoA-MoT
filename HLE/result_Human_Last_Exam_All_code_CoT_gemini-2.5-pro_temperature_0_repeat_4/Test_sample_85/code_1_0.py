import unicodedata

def reconstruct_sanskrit_ancestor():
    """
    This function reconstructs the Sanskrit ancestor for Hindi kãːs and Punjabi kə̃d͡ʒʱ.
    
    The reconstruction process is as follows:
    1.  k- vs. k-: The initial consonant is stable, so the ancestor starts with 'k'.
    2.  -ãːs vs. -ə̃d͡ʒʱ: The nasalized vowels point to a Vowel + Nasal sequence.
        The divergent consonants 's' and 'd͡ʒʱ' point to the Sanskrit cluster 'nṣ'.
        - Skt. 'anṣa' > Hindi 'ãːs' (simplification of ṣ > s, with compensatory lengthening and nasalization from n).
        - Skt. 'anṣa' > Punjabi 'ə̃d͡ʒʱ' (a regular development of the 'nṣ' cluster).
    3.  The most likely Sanskrit form is 'kanṣa'.
    4.  In IPA, this is transcribed as 'kɐnʂɐ'.
    """
    
    # The reconstructed Sanskrit word in Devanagari would be कंस or कन्स
    # The transliteration is 'kanṣa'
    # The IPA representation is 'kɐnʂɐ'
    
    reconstructed_ipa = "kɐnʂɐ"
    
    print("Reconstructed Sanskrit Ancestor (IPA):")
    print(reconstructed_ipa)

reconstruct_sanskrit_ancestor()