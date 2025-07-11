# Plan:
# 1. Identify the notes in the bass clef (bottom staff).
# 2. Assign beat values to each note type (quarter note = 1, half note = 2).
# 3. Sum the beat values to find the total beats in the measure.
# 4. Conclude the time signature from the total.

# Note durations in terms of quarter note beats.
quarter_note = 1
half_note = 2

# Notes found in the bass clef of the measure.
note_1_duration = quarter_note
note_2_duration = half_note
note_3_duration = quarter_note

# Calculate the total number of beats.
total_beats = note_1_duration + note_2_duration + note_3_duration

# Print the analysis and result.
print("To find the time signature, we sum the durations of the notes in the measure.")
print("The bass clef contains a quarter note, a half note, and another quarter note.")
print("The calculation for the total beats is:")
# The final equation with each number printed
print(f"{note_1_duration} + {note_2_duration} + {note_3_duration} = {total_beats}")
print(f"The measure contains a total of {total_beats} quarter-note beats.")
print("Therefore, the time signature is 4/4.")
