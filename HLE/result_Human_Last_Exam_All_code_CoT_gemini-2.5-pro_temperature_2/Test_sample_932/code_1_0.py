# The indices of tribes whose immatures are unlikely to be collected
# via the beat-sheet method due to their hidden or internal lifestyles.
unlikely_tribe_indices = [1, 3, 4, 7]

# The list is already in ascending order.
# We will convert the numbers to strings to print them separated by commas.
# This fulfills the request to output each number in the final answer.
output_string = ", ".join(map(str, unlikely_tribe_indices))

print(output_string)