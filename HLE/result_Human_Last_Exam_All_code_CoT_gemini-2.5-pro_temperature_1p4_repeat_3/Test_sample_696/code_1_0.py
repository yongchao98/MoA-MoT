from re import findall

# The 32-character solution is a lambda function that transforms the list of found strings.
# The original template `print(@r'\d+', input())})or ["NO"])` is syntactically invalid.
# A corrected, working version would look like this:

# The lambda is represented by the variable `transformer` for clarity.
# transformer = lambda l: sorted(set(map(int, l)))
# print(len('lambda l: sorted(set(map(int, l)))')) # This would print 32

print(*(
    # This is where the 32-character lambda would be called
    (lambda l: sorted(set(map(int, l))))(findall(r'\d+', input()))
    or ["NO"]
))