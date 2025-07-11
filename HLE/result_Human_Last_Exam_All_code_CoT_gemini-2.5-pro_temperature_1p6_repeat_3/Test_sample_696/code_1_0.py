import re

# The problem asks for a 32-character string to make a given one-liner work.
# The target one-liner, after fixing an apparent typo in the prompt's template, is:
# from re import*;print(* sorted({int(n)for n in findall(r'\d+', input())}) or ["NO"])
# The 32-character replacement for the placeholder '@' is: * sorted({int(n)for n in findall(
# This script breaks down the logic of that one-liner.

# Step 1: Get the input string from the user.
input_string = input()

# Step 2: Use regex's findall to extract all number strings.
# For "hello 123 world 456 hello 123", this gives ['123', '456', '123']
numbers_as_strings = re.findall(r'\d+', input_string)

# Step 3: Check if the resulting list is empty.
if not numbers_as_strings:
    # If no numbers were found, print "NO" as required.
    print("NO")
else:
    # Step 4: Use a set comprehension to convert strings to integers and get unique values.
    # For ['123', '456', '123'], this gives the set {123, 456}.
    unique_numbers_as_int = {int(n) for n in numbers_as_strings}

    # Step 5: Sort the unique numbers in ascending order.
    # The sorted() function returns a sorted list: [123, 456].
    sorted_numbers = sorted(unique_numbers_as_int)

    # Step 6: Print the sorted numbers.
    # As per the problem description, the final output should be the numbers separated by spaces.
    # We can use the join method for clarity here.
    # The list [123, 456] is converted to ['123', '456'] and then joined to "123 456".
    output_string = ' '.join(map(str, sorted_numbers))
    print(output_string)
