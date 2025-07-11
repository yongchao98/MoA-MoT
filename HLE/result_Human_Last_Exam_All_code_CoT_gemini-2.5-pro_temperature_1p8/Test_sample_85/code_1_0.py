import unicodedata

def reconstruct_sanskrit_ancestor():
    """
    This function analyzes the Hindi and Punjabi descendants to reconstruct the Sanskrit ancestor.
    
    1. Hindi 'kãːs' and Punjabi 'kə̃d͡ʒʱ' both point to an initial 'k' and a vowel 'a'.
    2. The nasalization in both (̃) suggests a lost nasal consonant from the ancestor, likely 'n'.
       So the form is 'kan...'.
    3. The divergence is in the final consonants: Hindi '-s' vs. Punjabi '-dʒʱ'.
    4. The Hindi '-s' (with preceding vowel lengthening and nasalization) strongly points to a 
       Sanskrit cluster like '-nṣ-' (n + retroflex sibilant ṣ).
       The change is kanṣa -> kaṃsa (Prakrit) -> kãːs (Hindi). This is a regular sound change.
    5. The Punjabi '-dʒʱ' is most regularly derived from a cluster like '-dhy-'. 
       If the ancestor were 'kandhya', it would yield 'kā̃jh' in Hindi, which is incorrect.
    6. Therefore, we must assume the ancestor is 'kanṣa' and that a more complex, though
       phonetically plausible, series of changes occurred in the evolution to Punjabi
       (e.g., palatalization, voicing, and aspiration of the '-nṣ-' cluster).
    7. The reconstructed Sanskrit word is 'kanṣa'.
    8. In IPA, this is represented as k + a + n + ṣ + a.
    """
    
    reconstruction_ipa = "kanʂa"
    
    # We can represent this with standard characters and add a note, or use the Unicode for ṣ.
    # IPA for ṣ is U+0282 (LATIN SMALL LETTER S WITH HOOK)
    
    sanskrit_k = "k"
    sanskrit_a1 = "a"
    sanskrit_n = "n"
    sanskrit_sh = "\u0282" # This is the IPA character 'ʂ'
    sanskrit_a2 = "a"

    # For clear output, we print each sound component.
    print(f"Reconstructed Sanskrit Ancestor (in IPA):")
    print(f"Initial sound: {sanskrit_k}")
    print(f"Vowel: {sanskrit_a1}")
    print(f"Nasal consonant: {sanskrit_n}")
    print(f"Main consonant (retroflex sibilant): {sanskrit_sh}")
    print(f"Final vowel: {sanskrit_a2}")
    
    final_equation = f"{sanskrit_k}{sanskrit_a1}{sanskrit_n}{sanskrit_sh}{sanskrit_a2}"
    print(f"\nFinal Equation: {final_equation}")

    # The final answer as a single string is returned here.
    return final_equation

# Execute the function and capture the final answer.
final_answer = reconstruct_sanskrit_ancestor()

# The final answer in the required format.
print(f"\n<<<{final_answer}>>>")
