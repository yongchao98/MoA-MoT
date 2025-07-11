import collections

# Step 1: Define the intervals for the ascending scale (Maqam Zanjaran).
# This is composed of Jins Hijaz on the root and Jins Rast on the fifth.
# Jins Hijaz (root): 0.5, 1.5, 0.5
# Connecting interval: 1.0
# Jins Rast (fifth): 1.0, 0.75, 0.75
ascending_intervals = [0.5, 1.5, 0.5, 1.0, 1.0, 0.75, 0.75]

# Step 2: Define the intervals for the descending scale.
# The descent from the octave to the 4th degree uses a modified scale
# with Jins Nahawand on the 4th degree.
# The descending notes are 8 -> 7 -> 6 -> 5 -> 4.
# The Jins Nahawand on the 4th gives the relative pitches of notes 4, 5, 6, and 7.
# Its pattern is Whole(1.0), Half(0.5), Whole(1.0).
# Interval from 8th (octave) down to 7th: 1.0
# Interval from 7th down to 6th (based on Jins Nahawand): 1.0
# Interval from 6th down to 5th (based on Jins Nahawand): 0.5
# Interval from 5th down to 4th (based on Jins Nahawand): 1.0
descending_intervals = [1.0, 1.0, 0.5, 1.0]

# Step 3: Combine the lists to get the full sequence of 11 intervals.
all_intervals = ascending_intervals + descending_intervals

# Step 4: Format and print the final output as requested.
# The format is {i1,i2,i3,...}
output_string = "{" + ",".join(map(str, all_intervals)) + "}"
print(output_string)

# The final answer in the required format.
final_answer = output_string
print(f"<<<{final_answer}>>>")