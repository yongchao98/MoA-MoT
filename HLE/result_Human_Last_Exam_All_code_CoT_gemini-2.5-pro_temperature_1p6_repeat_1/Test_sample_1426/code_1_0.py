import sys

def get_pitch_value(note_str):
    """
    Calculates a numerical pitch value for a given note string (e.g., 'F#5').
    This allows for easy comparison of pitches.
    """
    note_order = {
        'C': 0, 'B#': 0,
        'C#': 1, 'DB': 1,
        'D': 2,
        'D#': 3, 'EB': 3,
        'E': 4, 'FB': 4,
        'F': 5, 'E#': 5,
        'F#': 6, 'GB': 6,
        'G': 7,
        'G#': 8, 'AB': 8,
        'A': 9,
        'A#': 10, 'BB': 10,
        'B': 11, 'CB': 11
    }

    octave = int(note_str[-1])
    note_name = note_str[:-1].upper()

    if note_name not in note_order:
        raise ValueError(f"Unknown note name: {note_name}")

    base_value = note_order[note_name]
    return (octave * 12) + base_value

# Step 1: Identify the final three notes from the score.
# The last three notes for the first violin in Dvořák's Symphony No. 2 are
# in the final two measures of the 4th movement. They are F#5, A5, and D5.
final_notes = ['F#5', 'A5', 'D5']

print("The final three notes for the first violin are:")
# The prompt requires printing each "number" in the final equation.
# Here we print each note involved in the comparison.
for note in final_notes:
    print(note)

# Step 2: Find the lowest note among them.
lowest_note = None
min_pitch_value = sys.maxsize

for note in final_notes:
    pitch_value = get_pitch_value(note)
    if pitch_value < min_pitch_value:
        min_pitch_value = pitch_value
        lowest_note = note

print("\nCalculating the lowest note...")
print(f"The lowest of these notes is: {lowest_note}")
<<<D>>>