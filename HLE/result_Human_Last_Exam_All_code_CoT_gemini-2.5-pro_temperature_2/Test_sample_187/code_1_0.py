import json

def find_error_lines():
    # Line 11: The `^` operator on Python integers is bitwise XOR, not exponentiation.
    # This leads to an incorrect distance calculation.
    # Lines 22, 23: The `/` operator on Python integers performs float division in Python 3.
    # This violates the specification for "precise (non-floating point) arithmetic"
    # and will not produce the required exact rational results.
    error_lines = [11, 22, 23]
    error_lines.sort()
    print(error_lines)

find_error_lines()