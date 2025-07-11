# Define note durations as fractions of a whole note.
half_note = 0.5
quarter_note = 0.25

# The bottom staff contains two half notes.
note1_duration = half_note
note2_duration = half_note

# Calculate the total duration in the measure.
total_duration_in_whole_notes = note1_duration + note2_duration

# Determine the time signature based on the total duration.
# A total duration of 1.0 (a whole note) is equivalent to 4 quarter notes.
# This means there are 4 beats in the measure, and a quarter note gets one beat.
numerator = int(total_duration_in_whole_notes / quarter_note)
denominator = 4

print("Step 1: Identify the notes in the bottom staff.")
print("The measure contains one half note chord and one half note.")
print("Step 2: Assign numerical values to the note durations (where a whole note = 1).")
print(f"Value of a half note = {half_note}")
print("Step 3: Sum the durations of the notes in the measure.")
print(f"Total duration = {note1_duration} + {note2_duration} = {total_duration_in_whole_notes}")
print("Step 4: Interpret the result.")
print("A total duration of 1.0 is equivalent to one whole note, which is equal to four quarter notes.")
print(f"This means the time signature is {numerator}/{denominator}.")
