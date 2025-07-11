import sys
import io

# This script explains and generates the phonetic transcription
# for "shock absorber" as it would sound in Puerto Rican Spanish.

def generate_transcription():
    """
    Prints a step-by-step explanation and the final phonetic transcription.
    """
    print("This is a phonetic transcription of 'shock absorber' based on the rules of Puerto Rican Spanish phonology.")
    print("-" * 70)

    # --- Part 1: "Shock" ---
    print("Part 1: Adapting 'Shock' (/ʃɑk/)")
    print("1. /ʃ/ -> /tʃ/: The English 'sh' sound becomes the Spanish 'ch' sound.")
    print("2. /ɑ/ -> /o/: The English vowel sound is mapped to the Spanish /o/.")
    print("3. /k/ -> /k/: The final consonant sound remains the same.")
    shock_transcription = "[tʃok]"
    print(f"--> Result for 'Shock': {shock_transcription}\n")

    # --- Part 2: "Absorber" ---
    print("Part 2: Adapting 'Absorber' (/əbˈsɔːrbər/)")
    print("1. /bs/ -> /s/: The 'bs' cluster is simplified to 's', a common rule in Spanish (e.g., 'obscuro' -> 'oscuro'). The word is thus treated as 'asorber'.")
    print("2. /r/ -> /l/: The syllable-final /r/ is lateralized to an /l/ sound, a distinct feature of Puerto Rican Spanish. This happens in both '-sor-' and '-ber'.")
    print("3. /b/ -> [β]: The 'b' between two vowels is softened to a fricative sound.")
    absorber_transcription = "[a.solˈβel]"
    print(f"--> Result for 'Absorber': {absorber_transcription}\n")

    # --- Final Combination ---
    print("Part 3: Final Combination")
    print("The two adapted words are put together.")
    final_transcription = "[ˈtʃok a.solˈβel]"
    print("\nFinal Phonetic Equation:")
    print(f"{shock_transcription} + {absorber_transcription} = {final_transcription}")

generate_transcription()