# The problem of finding a path that visits every red ball exactly once is the Hamiltonian Path Problem.
# Let n be the number of red balls.

# Question A asks for the complexity of deciding if such a path exists.
# This is the decision version of the Hamiltonian Path problem.
# It is NP-complete, and the best-known algorithms have exponential complexity.
# A common complexity bound is O(2^n * n^2), which we simplify to O(2^n).
complexity_A = "O(2^n)"

# Question B asks for the complexity of finding the path, given it exists.
# This is the search version of the problem. It is NP-hard.
# The complexity of finding the path is the same as deciding its existence.
complexity_B = "O(2^n)"

# The final answer format is two complexities separated by a semicolon.
print(f"{complexity_A}; {complexity_B}")