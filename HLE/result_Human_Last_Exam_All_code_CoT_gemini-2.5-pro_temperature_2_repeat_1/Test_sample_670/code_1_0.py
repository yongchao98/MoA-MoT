# This script prints the Tate-style program-expression for the given code snippet.
# The expression represents the value of 'a' after the loop.
# It is constructed using the following logic:
# - The loop state is represented by a pair of variables (a, b),
#   corresponding to the program variables 'a' and 'i'.
# - The μ operator models the loop (cycle).
# - The λ operator defines the state update function within the loop.
# - The initial state for the loop computation is (1, 1).
# - The update function for state (a, b) is (a * b, b + 1).
# - The π₁ operator projects the final state pair to get the value of 'a'.

final_expression = "π₁(μ(λ(a, b).(a * b, b + 1)).(1, 1))"

print("The Tate-style program-expression is:")
print(final_expression)

# The numbers in the final equation are 1 (from π₁), 1 (from b + 1), 
# and 1, 1 from the initial state vector (1, 1).
print("\nThe numbers in the final equation are:")
print(1)
print(1)
print(1)
print(1)