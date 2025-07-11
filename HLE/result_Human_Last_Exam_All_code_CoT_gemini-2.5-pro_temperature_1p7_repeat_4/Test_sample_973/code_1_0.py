# The musical problem is to find the final melody note of "Happy Birthday" given a jazz chord progression.
# The melody of "Happy Birthday" ends on the tonic note. We need to find the final tonic.

# 1. Identify the final harmonic movement.
# The progression provided ends with Cm7 -> F7(9). This is a ii-V progression.
final_ii_chord = "Cm7"
final_v_chord = "F7(9)"
v_chord_root = "F"

# 2. A ii-V-I progression resolves to a tonic 'I' chord.
# The V chord's root is the 5th degree of the target key.
# We can find the tonic by moving up a perfect fourth (5 semitones) from the root of the V chord.

# Representing notes numerically for the calculation, where C=0, C#=1, D=2, etc.
notes = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B']
v_note_index = notes.index(v_chord_root)

# 3. Calculate the tonic note based on the V chord.
interval_to_tonic = 5  # A perfect fourth is 5 semitones
tonic_note_index = (v_note_index + interval_to_tonic) % 12
final_note = notes[tonic_note_index]

print(f"The final chord pair is {final_ii_chord} -> {final_v_chord}.")
print(f"This is a ii-V progression that resolves to an implied tonic (I).")
print(f"The root of the V chord is {v_chord_root}.")
print("To find the tonic, we ascend by a perfect fourth (5 semitones).")
print("\nHere is the final equation based on note indices (where C=0):")
# Print the equation with numbers
print(f"Tonic Index = (V Chord Root Index + Interval) mod 12")
print(f"Tonic Index = ({v_note_index} + {interval_to_tonic}) mod 12")
print(f"Tonic Index = {(v_note_index + interval_to_tonic)} mod 12 = {tonic_note_index}")
# The result
print(f"\nThe note at index {tonic_note_index} is {final_note}.")
print(f"Since the melody ends on the tonic, the note sung on the final 'you' is {final_note}.")