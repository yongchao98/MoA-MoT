import sys
import io

# The task is to reconstruct a Sanskrit ancestor from two invented Indo-Aryan words.
# This requires applying principles of historical linguistics. The code will outline
# the reconstructive logic and print the final result in IPA.

# Set output to UTF-8 to ensure IPA characters print correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def reconstruct_sanskrit_ancestor():
    """
    Performs and explains the reconstruction of a Sanskrit word from its
    Hindi and Punjabi descendants.
    """
    hindi_word = "kãːs"
    punjabi_word = "kə̃d͡ʒʱ"

    print("Reconstruction of the Sanskrit Ancestor")
    print("-" * 40)
    print(f"1. Descendant Words:")
    print(f"   - Hindi:   {hindi_word}")
    print(f"   - Punjabi: {punjabi_word}\n")

    print("2. Phonetic Analysis and Reconstruction:")

    # Step A: Initial Consonant
    initial_sanskrit = "k"
    print(f"   - Initial Consonant: Both words begin with 'k'. This is a very stable sound.")
    print(f"     Reconstructed initial: {initial_sanskrit}\n")

    # Step B: Vowel and Nasalization
    vowel_sanskrit = "a"
    print(f"   - Vowel: Hindi 'ãː' (long, nasalized) corresponding to Punjabi 'ə̃' (short, nasalized)")
    print(f"     points to an original short vowel '{vowel_sanskrit}' followed by a consonant cluster containing a nasal.")
    print(f"     The cluster's simplification caused compensatory lengthening in Hindi (a > aː).\n")

    # Step C: Medial Consonant Cluster
    cluster_sanskrit = "śn"
    print(f"   - Medial Consonant Cluster: The core puzzle is Hindi 's' vs. Punjabi 'd͡ʒʱ'.")
    print(f"     This points to divergent reflexes from a single Sanskrit cluster.")
    print(f"     The most probable source is '{cluster_sanskrit}'.\n")

    print("3. Derivation Paths:")
    print(f"   From Sanskrit *kaśna:")
    print(f"   A) Path to Punjabi '{punjabi_word}':")
    print(f"      Skt. 'śn' → Pkt. 'ñh' → Punjabi 'd͡ʒʱ'")
    print(f"      The palatal cluster 'ñh' was strengthened to a voiced aspirated affricate 'd͡ʒʱ'.\n")

    print(f"   B) Path to Hindi '{hindi_word}':")
    print(f"      Skt. 'śn' → Pkt. 'ss' → Hindi 's'")
    print(f"      A secondary Prakrit reflex simplifies 'śn' to 'ss', which becomes 's' in Hindi.")
    print(f"      The vowel lengthens ('a'→'aː') and undergoes spontaneous nasalization ('aː'→'ãː').\n")


    print("4. Final Reconstruction:")
    part1 = initial_sanskrit
    part2 = vowel_sanskrit
    part3 = cluster_sanskrit
    part4 = vowel_sanskrit # Thematic vowel
    reconstructed_word = f"{part1}{part2}{part3}{part4}"
    
    # Using 'equation' format as hinted by the prompt
    print(f"   Assembled Proto-Word: {part1} + {part2} + {part3} + {part4} = {reconstructed_word}\n")
    
    ipa_reconstruction = "kəɕnə"
    print(f"The most probable Sanskrit ancestor in IPA is: {ipa_reconstruction}")


if __name__ == "__main__":
    reconstruct_sanskrit_ancestor()
