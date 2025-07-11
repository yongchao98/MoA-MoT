# Define note durations relative to a whole note (where a whole note = 1)
quarter_note = 1/4
eighth_note = 1/8
sixteenth_note = 1/16
thirty_second_note = 1/32

# --- Analysis of the Right Hand (Top Staff) ---

# Beat 1 consists of one eighth note and four 32nd notes.
beat_1_duration = eighth_note + (4 * thirty_second_note)

# Beat 2 consists of four 32nd notes and two 16th notes.
beat_2_duration = (4 * thirty_second_note) + (2 * sixteenth_note)

# Beat 3 consists of four 16th notes.
beat_3_duration = 4 * sixteenth_note

# Beat 4 consists of four 16th notes.
beat_4_duration = 4 * sixteenth_note

# Calculate the total duration of the measure by summing the beats.
total_duration = beat_1_duration + beat_2_duration + beat_3_duration + beat_4_duration

# The numerator of the time signature is the number of beat units in the measure.
# The denominator tells us the beat unit is a quarter note (4).
# A total_duration of 1.0 represents a full measure of 4/4 time.

print("Let's calculate the total duration of the measure to find the time signature.")
print("The measure can be broken down into four main beats.")
print(f"Beat 1: 1 eighth note (1/8) + 4 thirty-second notes (4/32) = {eighth_note:.3f} + {4*thirty_second_note:.3f} = {beat_1_duration:.3f}")
print(f"Beat 2: 4 thirty-second notes (4/32) + 2 sixteenth notes (2/16) = {4*thirty_second_note:.3f} + {2*sixteenth_note:.3f} = {beat_2_duration:.3f}")
print(f"Beat 3: 4 sixteenth notes (4/16) = {4*sixteenth_note:.3f} = {beat_3_duration:.3f}")
print(f"Beat 4: 4 sixteenth notes (4/16) = {4*sixteenth_note:.3f} = {beat_4_duration:.3f}")
print("-" * 20)
print("Total duration calculation (in terms of quarter notes):")
print(f"{beat_1_duration:.2f} + {beat_2_duration:.2f} + {beat_3_duration:.2f} + {beat_4_duration:.2f} = {total_duration:.2f}")
print("This corresponds to 1 quarter note + 1 quarter note + 1 quarter note + 1 quarter note.")
print("The sum is 4 quarter notes.")
print("Therefore, the time signature is 4/4.")