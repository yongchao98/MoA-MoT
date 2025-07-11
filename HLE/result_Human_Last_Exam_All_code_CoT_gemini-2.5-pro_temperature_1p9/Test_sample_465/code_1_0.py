# Based on the analysis of statements about chromatic and orbital chromatic roots,
# we have identified the true statements.

# A. Real orbital chromatic roots are bounded by the greatest real chromatic root. -> True
# B. Chromatic roots may not be real. -> True
# C. Real chromatic roots may take on negative values. -> False
# D. Real chromatic roots may take on non-integer values. -> True
# E. Chromatic polynomials may have roots between 0 and 1. -> False

# We collect the letters corresponding to the true statements.
true_statements = ['A', 'B', 'D']

# The required output format is a sorted string without delimiters.
# We sort the list of characters alphabetically.
sorted_statements = sorted(true_statements)

# We join the sorted characters into a single string.
answer = "".join(sorted_statements)

# Print the final answer.
print(answer)