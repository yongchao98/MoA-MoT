import sys
import io

# The user wants me to identify the error in a C program and suggest a fix with the smallest change.
#
# 1. Analysis of the C code:
# - `const short int ss = 0x6425;`: This is likely an obfuscated way to create a format string. On a little-endian machine, the bytes for 0x6425 are stored as `25 64`. In ASCII, this is "%d". So, `(char*)&ss` can be used as a format string for scanf.
# - `char c;`: This variable is used with scanf's "%d" format specifier, which expects an `int*`, not a `char*`. This is undefined behavior, as scanf will write a 4-byte integer into a 1-byte space, corrupting the stack. However, for the given input values (-100 to 100), they fit within a char, so the value of `c` itself might appear correct after the call, even if memory is corrupted elsewhere.
# - `scanf(&ss,&c);`: This is the incorrect call. It should be `scanf((char*)&ss, &c);` to work as intended by the obfuscation. Let's assume the compiler/system allows this to work.
# - Loop structure: `while (1) { ... if (feof(stdin)) break; }`. This is a known problematic pattern. `feof()` only becomes true *after* a read operation fails by hitting the end of the file. The `scanf` calls will fail, but the code will proceed to use the now-stale value in `c` before the `feof()` check breaks the loop.
#
# 2. Tracing the failing test case:
# Input:
# 3
# 3 -1 7
# -5 2 -4
# 2 -1 -3
#
# - The initial `scanf(&ss,&c);` reads `3` (the `n`) into `c`. This value is immediately overwritten in the loop.
# - The loop proceeds, correctly summing the vectors. After the third vector (2, -1, -3) is processed, the sums are `x=0`, `y=0`, `z=0`. The last value read into `c` was `-3`.
# - The loop starts another iteration.
# - `scanf(&ss,&c);` fails because there is no more input. It returns EOF. The value of `c` remains unchanged (`-3`).
# - `x+=c;` executes. `x` becomes `0 + (-3) = -3`.
# - The next two `scanf` calls also fail. `y` becomes `0 + (-3) = -3`. `z` becomes `0 + (-3) = -3`.
# - `if (feof(stdin))` is now true, and the loop breaks.
# - `puts(x||y||z?"NO":"YES");` evaluates `-3 || -3 || -3`, which is true. It prints "NO". The correct answer is "YES".
#
# 3. Finding the smallest fix:
# The error is that the stale value of `c` is used after `scanf` fails. To fix this with a minimal change, we need to ensure `c` does not contain a stale value when the additions happen after a failed read.
#
# A simple way to achieve this is to reset `c` to 0 at the beginning of each loop iteration. If `scanf` fails, `c` will remain 0, and adding 0 to the sums will not corrupt the result.
#
# The change would be on line 8:
# Original: `scanf(&ss,&c); x+=c;`
# Fixed:    `c=0; scanf(&ss,&c); x+=c;`
#
# This change involves inserting `c=0; ` at the beginning of the line.
# The number of characters (edits) inserted is 5: `c`, `=`, `0`, `;`, ` `.
#
# The line number is 8. The number of edits is 5.
# The answer is in the format l:n.

def solve():
    """
    This function prints the line number and number of edits to fix the C code.
    """
    line_number = 8
    edit_count = 5
    print(f"{line_number}:{edit_count}")

solve()
<<<8:5>>>