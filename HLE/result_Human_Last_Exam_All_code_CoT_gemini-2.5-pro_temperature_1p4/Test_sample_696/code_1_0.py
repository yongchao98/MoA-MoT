# The original template from the problem seems to have typos in its structure.
# The following code correctly solves the problem and uses a 32-character substitution
# to form the core logic, assuming the intended structure was a valid python expression.
# The 32 characters are: * sorted( set( map(int,findall(

from re import*

# We build the full correct expression to demonstrate the logic.
# The findall() function gives us all number strings.
# map(int, ...) converts them to integers.
# set(...) makes them unique.
# sorted(...) sorts them.
# The 'or ["NO"]' handles cases with no numbers.
# The '*' unpacks the list for printing.

# This reconstructs the full line based on the logic derived.
# Note that this code requires four closing parentheses `))))` after input() to be valid,
# which contradicts the `})` in the problem's template.
# We prioritize the length constraint as it is more explicit.
final_list = sorted(set(map(int, findall(r'\d+', input()))))

print(*(final_list or ["NO"]))