# The plan is to generate a string representing the Tate-style program-expression for the given loop.
# This expression uses the μ-operator to define a recursive function that models the loop's execution.

# 1. Identify the numerical constants from the source code.
#    - The initial value of the accumulator 'a' is 1.
#    - The initial value of the loop counter 'i' is 1.
#    - The counter 'i' is incremented by 1 in each step.
initial_accumulator_value = 1
initial_counter_value = 1
counter_increment = 1

# 2. Define the bound variable names according to the user's instruction (`a` for the first, `b` for the second, etc.).
#    - The μ-operator binds the recursive function name. This is the first bound variable: 'a'.
#    - The λ-operator binds the function's arguments, which represent the loop's state.
#      - The first argument represents the accumulator. This is the second bound variable: 'b'.
#      - The second argument represents the loop counter. This is the third bound variable: 'c'.
mu_variable = "a"
accumulator_arg = "b"
counter_arg = "c"

# 3. Assemble the parts of the expression into a single string.
#    - The expression is a call to a recursive function defined by μ.
#    - The function takes the initial state (initial accumulator, initial counter) as arguments.
#    - The function body contains an 'if' statement for the loop condition (`c < n`).
#    - The 'then' branch is the recursive call with the updated state (`a(b * c, c + 1)`).
#    - The 'else' branch is the termination value (`b`).
expression_string = (
    f"(μ{mu_variable}."
    f"(λ{accumulator_arg}, {counter_arg}. "
    f"if {counter_arg} < n then "
    f"{mu_variable}({accumulator_arg} * {counter_arg}, {counter_arg} + {counter_increment}) "
    f"else {accumulator_arg}))"
    f"({initial_accumulator_value}, {initial_counter_value})"
)

# 4. Print the final constructed expression.
print(expression_string)