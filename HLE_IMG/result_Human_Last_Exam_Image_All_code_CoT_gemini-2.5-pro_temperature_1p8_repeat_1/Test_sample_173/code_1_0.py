# Let's represent the duration of musical notes as fractions.
# A whole note is 1, a half note is 1/2, a quarter note is 1/4, etc.

# The bass clef contains three notes.
# Note 1: A quarter note
note1_duration = 1/4
# Note 2: A quarter note
note2_duration = 1/4
# Note 3: A quarter note
note3_duration = 1/4

# Calculate the total duration by adding the durations of all notes in the measure.
total_duration = note1_duration + note2_duration + note3_duration

# The numerator of the time signature represents the number of beats per measure.
# The denominator represents the note value that gets one beat (4 for a quarter note).
# The total duration of 3/4 corresponds to a 3/4 time signature.
numerator = 3
denominator = 4

print("Let's analyze the bass clef to find the time signature.")
print(f"Note 1 is a quarter note, its duration is {int(1/note1_duration)}/{int(1/note1_duration)} of a beat, or {note1_duration} of a whole note.")
print(f"Note 2 is a quarter note, its duration is {int(1/note2_duration)}/{int(1/note2_duration)} of a beat, or {note2_duration} of a whole note.")
print(f"Note 3 is a quarter note, its duration is {int(1/note3_duration)}/{int(1/note3_duration)} of a beat, or {note3_duration} of a whole note.")
print("\nTo find the time signature, we sum the durations:")
print(f"Total duration = {note1_duration} + {note2_duration} + {note3_duration} = {total_duration}")
print(f"A total duration of {total_duration} means the measure contains the equivalent of {numerator} quarter notes.")
print(f"Therefore, the time signature is {numerator}/{denominator}.")
