# The initial value of the variable 'a' is 1.
initial_a = 1

# The update operation inside the loop is `a = a * i`.
# We represent this with a lambda function `λa. *(a, i)`, where 'a' is the
# bound variable representing the value from the previous iteration.
# The `μ` operator represents the fixed-point semantics of the loop.
# The entire loop's effect on 'a' is captured by `μ(λa. *(a, i))`.
# This is then applied to the initial value of 'a'.
# The variable `initial_a` is used here to explicitly show where the number in the equation comes from.
final_expression = f"(μ(λa. *(a, i)))({initial_a})"

print(final_expression)