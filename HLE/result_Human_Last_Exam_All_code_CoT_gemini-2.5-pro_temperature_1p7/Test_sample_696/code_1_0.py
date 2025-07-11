import re

# The user's input string
# For example: "hello 123 world 456 hello 123"
# Or: "no numbers here"
s = input()

# The one-liner in question is:
# from re import*;print(@r'\d+', input())})or ["NO"])

# The 32 characters to replace @ are: *(sorted(map(int, set(findall(  
# Let's break down the resulting code:
# print(*(sorted(map(int, set(re.findall(  r'\d+', s))))) or ["NO"]))

# 1. re.findall(r'\d+', s)
#    Finds all substrings consisting of one or more digits.
#    e.g., for "hello 123 world 456 hello 123", it returns ['123', '456', '123']
#    e.g., for "no numbers here", it returns []
numbers_as_strings = re.findall(r'\d+', s)

# 2. set(numbers_as_strings)
#    Creates a set from the list to get unique elements.
#    e.g., {'123', '456'}
#    e.g., set()
unique_strings = set(numbers_as_strings)

# 3. map(int, unique_strings)
#    Converts each string number to an integer. This returns a map object.
#    e.g., map object yielding 123, 456
#    e.g., map object yielding nothing
numbers_as_ints = map(int, unique_strings)

# 4. sorted(numbers_as_ints)
#    Sorts the integers in ascending order. This returns a list.
#    e.g., [123, 456]
#    e.g., []
sorted_numbers = sorted(numbers_as_ints)

# 5. sorted_numbers or ["NO"]
#    If sorted_numbers is not empty (it's "truthy"), the expression is sorted_numbers.
#    If sorted_numbers is empty `[]` (it's "falsy"), the expression becomes ["NO"].
#    e.g., [123, 456]
#    e.g., ["NO"]
final_list = sorted_numbers or ["NO"]

# 6. print(*final_list)
#    The * unpacks the list.
#    print(*[123, 456]) is equivalent to print(123, 456), which prints "123 456".
#    print(*["NO"]) is equivalent to print("NO"), which prints "NO".
print(*final_list)
