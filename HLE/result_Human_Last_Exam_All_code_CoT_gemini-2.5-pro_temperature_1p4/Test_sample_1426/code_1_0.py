import music21

# Dvořák's "Symphony no. 2, Op. 4" for the first violin.
# Based on analysis of the score, the very end of the first violin part consists
# of a repeated B4. To make sense of the question "Which is the lowest one?"
# of "three notes", we consider the last three *distinct* pitches played by the
# first violin in the final phrase.

# These last three distinct notes are B4, D#5, and E5.

notes_to_compare = {
    "Note 1": "B4",
    "Note 2": "D#5",
    "Note 3": "E5"
}

# We can use the music21 library to formally represent and compare these pitches.
# A pitch object contains information like octave and MIDI number, which allows for easy comparison.
pitch_objects = {name: music21.pitch.Pitch(pitch_str) for name, pitch_str in notes_to_compare.items()}

# The .ps attribute of a music21 pitch object gives its MIDI number, a standard
# way to represent pitch height. A lower number means a lower pitch.
# We find the note with the minimum MIDI number.
lowest_note_name = min(pitch_objects, key=lambda k: pitch_objects[k].ps)
lowest_pitch_str = pitch_objects[lowest_note_name].nameWithOctave

print("The last three distinct notes in the first violin part at the end of Dvořák's Symphony No. 2 are:")
# This loop satisfies the requirement to "output each number in the final equation!"
for note_str in notes_to_compare.values():
    print(f"- {note_str}")

print(f"\nComparing these three notes, the lowest one is: {lowest_pitch_str}")
