# In music notation, we can assign numerical values to note durations.
# Let's define the beat value of common notes, assuming a quarter note gets one beat.
quarter_note = 1
half_note = 2

# The bottom staff of the measure contains one quarter note and one half note.
# We sum their values to find the total number of beats in the measure.
notes_in_measure = [quarter_note, half_note]
total_beats = sum(notes_in_measure)

# The total number of beats gives us the numerator (top number) of the time signature.
numerator = total_beats

# Since we are using the quarter note as our unit of one beat, the denominator is 4.
denominator = 4

# Print the calculation and the result.
print(f"Calculation: {notes_in_measure[0]} (quarter note) + {notes_in_measure[1]} (half note) = {total_beats} beats")
print(f"The time signature is {numerator}/{denominator}.")
