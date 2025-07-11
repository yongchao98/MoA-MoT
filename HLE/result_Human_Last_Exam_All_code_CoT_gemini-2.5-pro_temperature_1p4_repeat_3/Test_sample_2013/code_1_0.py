# Define the known information about the piece.
composer_surname = "Khachaturian"
opus_number = 38
start_measure = 813

# The duration of the 'off-beat' period is given in the problem description.
off_beat_duration_in_measures = 2

# Calculate the end measure of the musical joke.
# The joke starts in the start_measure, continues for the duration,
# and the recovery happens in the measure immediately following.
# However, the prompt asks for a range from the beginning of the flourish
# to the measure where the recovery occurs. Based on a 3-measure event
# structure (flourish/error, off-beat, recovery), the end measure is
# start_measure + 2.
end_measure = start_measure + off_beat_duration_in_measures

# Format the range as a hyphenated string.
measure_range = f"{start_measure}-{end_measure}"

# Print the final answer in the specified comma-separated format.
print(f"{composer_surname}, {opus_number}, {measure_range}")

# As requested, output the numbers in the equation used to find the end measure.
print("Equation for the end measure:")
print(f"{start_measure} + {off_beat_duration_in_measures} = {end_measure}")