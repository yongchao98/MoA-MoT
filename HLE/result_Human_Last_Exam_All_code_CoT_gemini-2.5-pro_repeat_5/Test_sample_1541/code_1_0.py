# Read the input sum as a string.
# Example input: "3+2+1"
s = input()

# The summands are separated by '+'. We can split the string by '+'
# to get a list of the numbers (as strings).
# For "3+2+1", this results in ['3', '2', '1']
summands_as_strings = s.split('+')

# To sort the numbers correctly, we first need to convert them from strings to integers.
numbers = []
for summand in summands_as_strings:
    numbers.append(int(summand))

# Now, we sort the list of integers.
# For [3, 2, 1], this becomes [1, 2, 3]
numbers.sort()

# To print the final sum in the required format (e.g., "1+2+3"),
# we need to convert the sorted numbers back into strings.
output_list = []
for number in numbers:
    output_list.append(str(number))

# Finally, we join the elements of our list of strings with a '+' sign in between.
# For ['1', '2', '3'], this becomes "1+2+3"
final_sum_string = "+".join(output_list)

# Print the rearranged sum that Xenia can calculate.
print(final_sum_string)