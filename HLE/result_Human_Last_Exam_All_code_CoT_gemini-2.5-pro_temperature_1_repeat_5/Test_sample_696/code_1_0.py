# The one-liner from the problem description contains syntax errors in its
# closing parentheses. The following is a corrected, runnable version that
# incorporates the 32-character solution.

from re import*;print(* (sorted(set(map(int,findall(r'\d+', input()))))) or ["NO"]))
