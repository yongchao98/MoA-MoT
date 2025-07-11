from re import *

# The user's provided snippet:
# print(@r'\d+', input())})or ["NO"])
# where '@' is the 32-character string we found.

# We will construct and execute the full line of code to solve the problem.
# The core logic is to find all numbers, convert them to integers for sorting,
# make them unique with a set, sort them, and then print them space-separated.
# If no numbers are found, "NO" is printed.

try:
    # We find all sequences of digits, which gives a list of strings.
    numbers = findall(r'\d+', input())

    # We process the list of strings if it's not empty.
    if numbers:
        # 1. Convert strings to integers: `int(n) for n in numbers`
        # 2. Use a set for uniqueness: `{...}`
        # 3. Sort numerically: `sorted({...})`
        unique_sorted_numbers = sorted({int(n) for n in numbers})

        # The problem asks for space-separated output.
        # We can convert each number back to a string and join.
        print(' '.join(map(str, unique_sorted_numbers)))
    else:
        # If no numbers were found, print "NO".
        print("NO")

except EOFError:
    # This is to handle cases where input is not provided,
    # e.g., in an automated testing environment without redirection.
    print("NO")
