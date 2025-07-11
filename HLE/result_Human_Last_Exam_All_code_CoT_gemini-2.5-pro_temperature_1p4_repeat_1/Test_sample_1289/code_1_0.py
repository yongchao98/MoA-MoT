# 1. Define the intervals for the ascending Maqam Zanjaran.
# The musician sings 8 notes, creating 7 intervals.
# Jins Ajam on tonic: Whole (1), Whole (1), Semi (0.5)
# Connecting tone: Whole (1)
# Jins Hijaz on fifth: Semi (0.5), Augmented Second (1.5), Semi (0.5)
ascending_intervals = [1, 1, 0.5, 1, 0.5, 1.5, 0.5]

# 2. Define the intervals for the descending modified scale.
# The musician descends from note 8 to note 4, creating 4 intervals.
# The modification is Jins Nahawand on the 4th degree.
# Jins Nahawand: Whole (1), Semi (0.5), Whole (1). This defines the scale from note 4 to 7.
# An additional whole tone (1) connects the 7th to the octave.
# The upper part of the scale is thus: Note 4 ->(1)-> N5 ->(0.5)-> N6 ->(1)-> N7 ->(1)-> Octave(N8)
# The musician descends from N8 to N4: N8 -> N7 -> N6 -> N5 -> N4.
# The intervals are:
# Interval from N8 down to N7: 1
# Interval from N7 down to N6: 1
# Interval from N6 down to N5: 0.5
# Interval from N5 down to N4: 1
descending_intervals = [1, 1, 0.5, 1]

# 3. Combine the lists to get all 11 intervals in order.
all_intervals = ascending_intervals + descending_intervals

# 4. Format the result as a string and print it.
# Convert each number to its string representation to handle floats correctly.
# The problem asks to output each number in the final equation.
# We will construct the string as requested.
final_output_string = "{" + ",".join(map(str, all_intervals)) + "}"

print(final_output_string)
<<<{"1,1,0.5,1,0.5,1.5,0.5,1,1,0.5,1"}>>>