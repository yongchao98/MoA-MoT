import unicodedata

def get_ipa_char_names(ipa_string):
    """Prints the name of each character in the IPA string."""
    print(f"Reconstructed Sanskrit ancestor in IPA: {ipa_string}")
    print("\nCharacter breakdown:")
    for char in ipa_string:
        try:
            name = unicodedata.name(char)
            print(f"'{char}': {name}")
        except ValueError:
            print(f"'{char}': Basic ASCII character")

# The reconstructed Sanskrit ancestor for Hindi kãːs and Punjabi kə̃d͡ʒʱ is kaṁkṣa.
# We will now print its representation in the International Phonetic Alphabet (IPA).
# k = k
# a = a (representing the short 'a' vowel, phonetically [ɐ])
# ṁ (anusvara before a velar stop) = ŋ (voiceless velar nasal)
# kṣ (cluster) = kʂ (k + voiceless retroflex fricative)
# a = a
reconstructed_ipa = "kaŋkʂa"

get_ipa_char_names(reconstructed_ipa)
