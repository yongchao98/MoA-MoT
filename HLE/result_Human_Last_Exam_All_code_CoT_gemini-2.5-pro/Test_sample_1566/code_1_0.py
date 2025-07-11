# This script calculates the number of topologically distinct continua
# based on the logical deduction outlined above.

# Case 1: The continuum is composed of end points (E) and interior points (I).
# This structure corresponds to a topological arc (homeomorphic to [0,1]).
# We verify if it meets the given properties.
# Property (1): An arc has 2 end points. Since 1 < 2 < infinity, this is satisfied.
# Property (2): An arc has exactly two orbits under its homeomorphism group:
#   - Orbit 1: The set of its two endpoints, {0, 1}.
#   - Orbit 2: The set of its interior points, (0, 1).
# This condition is satisfied.
# Therefore, the arc is a valid solution.
solutions_from_case_arc = 1

# Case 2: The continuum is composed of end points (E) and branch points (B).
# A space with only end points and branch points cannot be connected
# without interior points to form edges. Thus, it cannot be a continuum.
solutions_from_case_branching = 0

# The total number of topologically distinct continua is the sum of solutions
# from all possible valid cases.
total_solutions = solutions_from_case_arc + solutions_from_case_branching

print("The total number of topologically distinct continua satisfying the properties is calculated as follows:")
print("Number of solutions from arc-like structures without branch points = 1")
print("Number of solutions from structures with branch points but no interior points = 0")
print("\nThe final equation is the sum of the number of solutions from each valid case:")
# Outputting each number in the final equation as requested.
print(solutions_from_case_arc, "+", solutions_from_case_branching, "=", total_solutions)

# The final answer is the result of this calculation.
print(f"\nThus, there is only {total_solutions} such continuum.")