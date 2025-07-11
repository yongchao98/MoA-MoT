# The final code is constructed by replacing @ with the 32 characters.
# Although the original template has a minor bracket mismatch,
# this code implements the intended logic to solve the problem.
#
# Placeholder (@): *sorted({int(i) for i in findall 
#
# Assembling the line:
# print(*sorted({int(i) for i in findall (r'\d+', input())}) or ["NO"])
# Note: For the code to be syntactically correct, an extra ')' is needed at the end.
# print(*(sorted({int(i) for i in findall(r'\d+', input())}) or ["NO"]))

# Here is the code to execute:
from re import *

# The problem is to find the 32-character string for the placeholder.
# The code below is the fully formed, correct one-liner.
# It finds all number strings, converts them to a set of unique integers,
# sorts them, and then uses the splat (*) operator to print them space-separated.
# The `or ["NO"]` trick handles the case where no numbers are found.
print(*(sorted({int(i) for i in findall(r'\d+', input())}) or ["NO"]))