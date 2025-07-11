# Define the duration of a sixteenth note as a fraction of a whole note.
sixteenth_note = 1/16

# Calculate the duration of each beat in the measure.
# Beat 1: four sixteenth notes
beat_1_duration = 4 * sixteenth_note
# Beat 2: one sixteenth rest and three sixteenth notes
beat_2_duration = 1 * sixteenth_note + 3 * sixteenth_note
# Beat 3: four sixteenth notes
beat_3_duration = 4 * sixteenth_note
# Beat 4: one sixteenth rest and three sixteenth notes
beat_4_duration = 1 * sixteenth_note + 3 * sixteenth_note

# Sum the durations of the beats to get the total duration of the measure.
total_duration = beat_1_duration + beat_2_duration + beat_3_duration + beat_4_duration

# Express the result in terms of quarter notes.
# A quarter note has a duration of 1/4.
quarter_note = 1/4
num_quarter_notes = total_duration / quarter_note

print("Analyzing the rhythm of the measure:")
print(f"Beat 1 consists of 4 sixteenth notes. Duration = 4 * 1/16 = {beat_1_duration:.2f}")
print(f"Beat 2 consists of 1 sixteenth rest + 3 sixteenth notes. Duration = 1/16 + 3/16 = {beat_2_duration:.2f}")
print(f"Beat 3 consists of 4 sixteenth notes. Duration = 4 * 1/16 = {beat_3_duration:.2f}")
print(f"Beat 4 consists of 1 sixteenth rest + 3 sixteenth notes. Duration = 1/16 + 3/16 = {beat_4_duration:.2f}")
print("\nCalculating the total duration:")
print(f"Total duration = {beat_1_duration:.2f} + {beat_2_duration:.2f} + {beat_3_duration:.2f} + {beat_4_duration:.2f} = {total_duration}")
print(f"A whole note has a duration of 1. A quarter note has a duration of {quarter_note:.2f}.")
print(f"The measure contains {int(num_quarter_notes)} quarter-note beats.")
print("\nThis means the time signature has a '4' on top (4 beats) and a '4' on the bottom (quarter note gets the beat).")
print("The time signature is 4/4.")