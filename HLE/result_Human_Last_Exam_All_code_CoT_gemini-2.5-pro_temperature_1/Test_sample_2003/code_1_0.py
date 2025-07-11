import sys

def solve_key_signature_formula():
    """
    This function derives the formula for the sum of sharps in key signatures
    based on the problem's rules.
    """

    # Step 1: Establish base values for natural notes.
    # The number of sharps in the major key signatures for C,D,E,G,A,B.
    base_sharps = {'C': 0, 'D': 2, 'E': 4, 'G': 1, 'A': 3, 'B': 5}
    
    # Per the problem, F major is treated as E# major.
    # S(E#) = S(E) + 7 sharps.
    sharps_F = base_sharps['E'] + 7
    
    # Step 2 & 3: Calculate the number of sharps for each of the 12 initial notes
    # and find their sum (the constant term T0).
    initial_notes = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    sharps_n0 = []
    
    for note in initial_notes:
        base_note = note[0]
        num_sharps_in_name = len(note) - 1

        if base_note == 'F':
            # For F: use the calculated S(F)=11. For F#: S(F#) = S(F) + 7.
            sharps = sharps_F + num_sharps_in_name * 7
        else:
            # For all other notes: S(Note#) = S(Note) + 7
            sharps = base_sharps[base_note] + num_sharps_in_name * 7
        sharps_n0.append(sharps)
        
    T0 = sum(sharps_n0)

    # Step 4: Calculate the coefficient for n.
    # Each of the 12 notes is sharped 'n' times. Each sharp adds 7 to the key signature.
    # Total addition = 12 * 7 * n
    n_coefficient = 12 * 7

    # Step 5: Print the derivation and final formula.
    print("Derivation of the Formula:")
    print("-" * 30)
    print(f"1. Base sharps for natural notes (with F as E#):")
    print(f"   C: {base_sharps['C']}, D: {base_sharps['D']}, E: {base_sharps['E']}, F(E#): {sharps_F}, G: {base_sharps['G']}, A: {base_sharps['A']}, B: {base_sharps['B']}\n")

    print(f"2. Number of sharps for all 12 initial notes (for n=0):")
    note_sharps_str = ", ".join([f"{note}: {s}" for note, s in zip(initial_notes, sharps_n0)])
    print(f"   [{note_sharps_str}]\n")

    print(f"3. The sum for n=0 (the constant term):")
    sharps_sum_str = " + ".join(map(str, sharps_n0))
    print(f"   Sum = {sharps_sum_str}")
    print(f"   Constant Term = {T0}\n")
    
    print(f"4. The coefficient for the 'n' term:")
    print(f"   Each of the 12 keys gets 7*n additional sharps.")
    print(f"   Total additional sharps = 12 * 7 * n")
    print(f"   'n' Coefficient = {n_coefficient}\n")
    
    print("5. The Final Formula S(n) = (n_coefficient * n) + (constant_term):")
    # This line prints the final equation with the calculated numbers as requested.
    print(f"\nThe final formula is: S(n) = {n_coefficient}n + {T0}")

solve_key_signature_formula()