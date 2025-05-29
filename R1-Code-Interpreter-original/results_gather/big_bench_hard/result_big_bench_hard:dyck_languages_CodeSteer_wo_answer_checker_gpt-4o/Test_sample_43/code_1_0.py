# Input sequence
sequence = "{ { { } }"

# Count the number of opening and closing braces
opening_braces = sequence.count('{')
closing_braces = sequence.count('}')

# Calculate the number of closing braces needed
needed_closing_braces = opening_braces - closing_braces

# Complete the sequence by adding the needed closing braces
completed_sequence = sequence + '}' * needed_closing_braces

# Output the final answer
print(f"<<<{completed_sequence}>>>")