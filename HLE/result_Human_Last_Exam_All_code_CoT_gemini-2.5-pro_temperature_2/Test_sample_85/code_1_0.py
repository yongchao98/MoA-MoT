import sys
# This script reconstructs the Sanskrit ancestor of two Indo-Aryan words based on phonological rules.
# Descendant forms:
# Hindi: kãːs
# Punjabi: kə̃d͡ʒʱ

# Based on linguistic analysis, the reconstructed components of the Sanskrit word in IPA are as follows:

# Initial consonant: 'k' is consistent in both descendants.
p1 = "k"
# Vowel: Sanskrit short 'a', phonetically a schwa.
p2 = "ə"
# Medial cluster: Reconstructed as 'msy' based on the s / d͡ʒʱ correspondence.
p3 = "m"
p4 = "s"
p5 = "j"
# Final vowel: Standard thematic vowel, short 'a'.
p6 = "ə"

# To satisfy the instruction "output each number in the final equation",
# we will print each phoneme of the reconstructed word.
print(f"Phoneme 1: {p1}")
print(f"Phoneme 2: {p2}")
print(f"Phoneme 3: {p3}")
print(f"Phoneme 4: {p4}")
print(f"Phoneme 5: {p5}")
print(f"Phoneme 6: {p6}")

# The full reconstructed word in IPA is the concatenation of these phonemes.
reconstructed_ipa = p1 + p2 + p3 + p4 + p5 + p6

# Print the final reconstructed string.
print("\nFinal Reconstructed String (IPA):")
print(reconstructed_ipa)

# Redirecting stderr to a dummy stream to avoid any potential messages from libraries if this were a more complex script.
# This ensures the final output for the user is clean.
class DevNull:
    def write(self, msg):
        pass
sys.stderr = DevNull()
