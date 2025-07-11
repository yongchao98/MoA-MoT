# This script prints the ordered list of line numbers that contain errors.
# The errors stem from the fact that an imported Python file is not processed by the Sage preparser,
# so standard Python rules for operators apply.
# '^' is XOR, not power, and '/' is float division, not exact rational division.

error_lines = [11, 18, 22, 23, 29, 30, 31]
print(error_lines)