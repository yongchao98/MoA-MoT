#
# Step 1: Analyze the properties of the formula φ.
#
# The formula φ has n distinct atomic variables (n >= 2).
#
# Condition 1: φ has 2^(n-1) satisfying truth-value assignments.
# Total possible assignments for n variables = 2^n.
# Satisfying assignments = 2^(n-1) = (2^n) / 2.
# This means the formula is true for exactly half of all possible inputs.
#
# Condition 2: φ is not a tautology.
# A tautology is true for all 2^n assignments. Since φ is true for only half, this
# condition is automatically satisfied by Condition 1.
#
# Condition 3: Any two distinct assignments v1 and v2 that make φ true must differ on
# at least one atomic variable.
# By definition, any two distinct assignments v1 and v2 must differ on at least one
# variable. This condition is redundant and provides no new constraints on φ.
#
# Step 2: Determine the minimum complexity for such a formula.
#
# The problem asks for the minimum number of variables in any logically equivalent
# formula ψ. This is equivalent to asking for the minimum number of variables that
# a function with these properties must depend on.
#
# Let's see how many variables such a function must depend on.
#
# Can it depend on 0 variables?
# A function of 0 variables is a constant.
# - If constant TRUE: It would be a tautology, violating Condition 2.
# - If constant FALSE: It would have 0 satisfying assignments. Since n>=2, 2^(n-1) is at least 2.
#   So it cannot be the constant FALSE function.
# Conclusion: The function must depend on at least one variable.
#
# Can it depend on just 1 variable?
# Let's consider a function of n variables that only depends on the first one, p_1.
# Example function f(p_1, p_2, ..., p_n) = p_1.
# - This function is true whenever p_1 is TRUE.
# - The other n-1 variables can be anything (2^(n-1) possibilities).
# - So, this function is true for exactly 2^(n-1) assignments. This satisfies Condition 1.
# - This function is not a tautology (it's false when p_1 is FALSE). This satisfies Condition 2.
#
# Step 3: Conclude the minimum number of variables.
#
# We have found that a function satisfying the given conditions can depend on as few as one variable.
# A formula φ could be constructed to be equivalent to this simple function (e.g., φ = p_1 ∧ (p_2 ∨ ¬p_2) ∧ ...).
# An equivalent formula ψ would represent the same function. The most simplified form of this function is just "p_1".
# This simplified formula, ψ = p_1, uses only one variable.
# We already established that the function must depend on at least one variable.
# Therefore, the minimum number of distinct atomic variables required in ψ is 1.
# The available connectives {∧, ¬} are functionally complete and can represent any function,
# including this simple one.

min_variables = 1

print("The logic demonstrates that a function satisfying the conditions can depend on as few as one variable.")
print("The minimal representation of such a function would therefore require only that one variable.")
print("The minimum number of distinct atomic variables required is:")
print(min_variables)
