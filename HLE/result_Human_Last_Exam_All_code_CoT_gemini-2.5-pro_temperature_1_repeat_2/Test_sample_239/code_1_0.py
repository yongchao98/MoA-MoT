import math

# Input values from the question
n = 4 * 10**9
m = 4 * 10**9
a = 1

# The number of flagstones for a single dimension of length x and flagstone size a
# is ceil(x/a). Using integer arithmetic, this is equivalent to (x + a - 1) // a.
# This method avoids floating-point inaccuracies and is generally more efficient.

flagstones_n = (n + a - 1) // a
flagstones_m = (m + a - 1) // a

total_flagstones = flagstones_n * flagstones_m

# The problem requires printing the final number of flagstones.
print(int(total_flagstones))