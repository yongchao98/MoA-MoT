# Plan:
# 1. Define the intervals for the ascending Maqam Zanjaran scale.
#    - Maqam Zanjaran consists of Jins Zanjaran on C and Jins Hijazkar on G.
#    - Jins Zanjaran (C to F): Intervals are {0.75, 1.25, 0.5}
#    - Linking interval (F to G): {1.0}
#    - Jins Hijazkar (G to C): Intervals are {0.5, 1.5, 0.5}
ascending_intervals = [0.75, 1.25, 0.5, 1.0, 0.5, 1.5, 0.5]

# 2. Define the intervals for the descending path.
#    - The descent is from the octave (C) to the 4th degree (F).
#    - The scale uses Jins Nahawand on the 4th degree (F).
#    - Jins Nahawand on F (F, G, Ab, Bb) has intervals {1.0, 0.5, 1.0}.
#    - The note before the octave C is Bb, making the interval from Bb to C a whole tone (1.0).
#    - The sung descending notes are C -> Bb -> Ab -> G -> F.
#    - Interval from C to Bb: 1.0
#    - Interval from Bb to Ab: 1.0
#    - Interval from Ab to G: 0.5
#    - Interval from G to F: 1.0
descending_intervals = [1.0, 1.0, 0.5, 1.0]

# 3. Combine the ascending and descending intervals.
all_intervals = ascending_intervals + descending_intervals

# 4. Format the list into the required string format "{n1,n2,n3...}".
#    - Convert each number in the list to its string representation.
#    - Join the string representations with a comma.
#    - Add the curly braces at the beginning and end.
output_string = "{" + ",".join(map(str, all_intervals)) + "}"

# 5. Print the final result.
print(output_string)

# The final list of intervals is {0.75, 1.25, 0.5, 1.0, 0.5, 1.5, 0.5, 1.0, 1.0, 0.5, 1.0}