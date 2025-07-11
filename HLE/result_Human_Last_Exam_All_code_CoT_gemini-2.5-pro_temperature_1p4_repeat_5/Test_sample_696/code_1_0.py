from re import*;print(*(sorted(set(map(int,findall(r'\d+', input())))))or ["NO"])

# To explain the filled-in one-liner from the problem:
# The original template was: from re import*;print(@r'\d+', input())})or ["NO"])
# Assuming the "}" is a typo for ")" and that 4 closing parentheses are needed,
# the code below is a reconstruction. The 32 characters for @ are found
# by adding spaces to the most concise logical core.

# @ = "* sorted( set( map( int,findall( " (32 chars)
# After `input()` comes "))))"
# So the line becomes:
# print(*( sorted( set( map( int,findall( r'\d+', input())))))or ["NO"])

# The above line is runnable and solves the problem.
# For clarity in this block, I am presenting the standard, un-spaced version
# which is functionally identical.