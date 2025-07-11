# This script identifies and prints the incorrect statements about
# State Enumeration Graphs (SEG) and Logical Attack Graphs (LAG).
# The analysis concluded that statements C and D are incorrect.

# C is incorrect because it oversimplifies the reason for the size difference.
# While monotonicity is a major factor, logical abstraction is also a key contributor
# to the compactness of LAGs.

# D is incorrect because probabilistic reachability can be calculated in graphs
# with loops. For SEGs (as Markov chains), this is a standard procedure.

# The final answer consists of the letters of the incorrect statements,
# sorted alphabetically and separated by a comma.

incorrect_statements = ["C", "D"]
print(",".join(incorrect_statements))