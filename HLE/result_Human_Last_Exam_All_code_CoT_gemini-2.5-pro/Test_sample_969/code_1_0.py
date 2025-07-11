# The pattern of the sequence is based on runs of numbers.
# By analyzing the structure, we can deduce the underlying rules for generating the sequence.

# The sequence of numbers being repeated is [3, 2, 1, 2, 3, 2, 1].
# The sequence of the lengths of these runs is [1, 1, 1, 1, 3, 3, 3].

# We can generate the full sequence from these two patterns.
# The given sequence is the first 9 elements of this fully generated sequence.
# We will generate the sequence and then print out the first 13 elements.

numbers_pattern = [3, 2, 1, 2, 3, 2, 1]
lengths_pattern = [1, 1, 1, 1, 3, 3, 3]

full_sequence = []
for i in range(len(numbers_pattern)):
    number = numbers_pattern[i]
    length = lengths_pattern[i]
    for _ in range(length):
        full_sequence.append(number)

# The problem asks to complete the sequence.
# The initial sequence given has 9 elements.
# The next 4 elements will be elements at index 9, 10, 11, and 12 of our generated sequence.
final_sequence_to_print = full_sequence[:13]

# Print each number in the completed sequence, as requested.
print(*final_sequence_to_print)