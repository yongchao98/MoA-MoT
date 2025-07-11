# Read the input sum as a string
s = input()

# Extract the numbers by splitting the string by the '+' delimiter.
# For example, "1+3+2+1" becomes ['1', '3', '2', '1']
summands = s.split('+')

# Sort the list of numbers. Because they are single digits,
# a standard lexicographical sort works perfectly.
# ['1', '3', '2', '1'] becomes ['1', '1', '2', '3']
summands.sort()

# Join the sorted list of numbers back into a string,
# separated by the '+' character.
# ['1', '1', '2', '3'] becomes "1+1+2+3"
new_sum = "+".join(summands)

# Print the final rearranged sum.
print(new_sum)