# Plan:
# 1. Define the intervals for the ascending scale (Maqam Zanjaran).
# 2. Define the intervals for the descending scale (modified with Jins Nahawand).
# 3. Combine the interval lists and print in the specified format.

# Step 1: Ascending Intervals
# Maqam Zanjaran on C is composed of Jins Rast on C and Jins Zanjaran on G.
# Scale notes: C, D-half-flat, E-half-flat, F, G, A-flat, B, C
# Interval calculation:
# C to D-half-flat: 0.75 (three-quarter tone)
# D-half-flat to E-half-flat: 0.75 (three-quarter tone)
# E-half-flat to F: 1.0 (whole tone)
# F to G: 1.0 (whole tone) - connecting interval
# G to A-flat: 0.5 (semitone)
# A-flat to B: 1.5 (one and a half tones)
# B to C: 0.5 (semitone)
ascending_intervals = [0.75, 0.75, 1.0, 1.0, 0.5, 1.5, 0.5]

# Step 2: Descending Intervals
# The descent uses a scale modified with Jins Nahawand on the 4th degree (F).
# This replaces the upper part of the scale, changing B natural to B-flat.
# Modified scale for descent: C, D-half-flat, E-half-flat, F, G, A-flat, B-flat, C.
# The musician descends from the octave (C) to the 4th degree (F).
# Notes sung: C -> B-flat -> A-flat -> G -> F.
# Interval calculation:
# C down to B-flat: 1.0 (whole tone)
# B-flat down to A-flat: 1.0 (whole tone)
# A-flat down to G: 0.5 (semitone)
# G down to F: 1.0 (whole tone)
descending_intervals = [1.0, 1.0, 0.5, 1.0]

# Step 3: Combine and Print
# The total performance has 11 intervals: 7 ascending + 4 descending.
total_intervals = ascending_intervals + descending_intervals

# Format the list as a string "{val1,val2,...}"
# The problem asks to output each number in the final equation.
# So we will format the string for the final print statement.
output_string = "{" + ",".join(map(str, total_intervals)) + "}"

print(output_string)
