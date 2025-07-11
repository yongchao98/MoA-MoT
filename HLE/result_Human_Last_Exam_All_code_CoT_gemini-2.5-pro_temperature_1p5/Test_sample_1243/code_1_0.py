import math

# (a) Calculate the number of non-zero binary vectors of length 16 with weight divisible by 4.
n_a = 16
# These are weights 4, 8, 12, 16.
c16_4 = math.comb(n_a, 4)
c16_8 = math.comb(n_a, 8)
c16_12 = math.comb(n_a, 12)
c16_16 = math.comb(n_a, 16)
ans_a = c16_4 + c16_8 + c16_12 + c16_16
eq_a = f"{c16_4} + {c16_8} + {c16_12} + {c16_16} = {ans_a}"

# (b) The answer is determined by a theoretical argument.
ans_b = "no"

# (c) Calculate the number of non-zero binary vectors of length 12 with weight divisible by 4.
n_c = 12
# These are weights 4, 8, 12.
c12_4 = math.comb(n_c, 4)
c12_8 = math.comb(n_c, 8)
c12_12 = math.comb(n_c, 12)
ans_c = c12_4 + c12_8 + c12_12
eq_c = f"{c12_4} + {c12_8} + {c12_12} = {ans_c}"

# Print the final answer in the required format.
# The user wants each number in the final equation.
# This print statement is for the final answer block.
final_answer_str = f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}"
print(f"The final answer is formatted as required.\n(a) is calculated from C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16)\n(c) is calculated from C(12, 4) + C(12, 8) + C(12, 12)")
# This print is not the final answer. The requested format is below
# Let's adjust to be very direct
# We will construct a string to be wrapped by the special <<<>>> tags later.
final_string = f"(a) {eq_a}; (b) {ans_b}; (c) {eq_c}"
# This will be used for the final output. The problem asks me to print though.
# And only one code block.
print(f"({c16_4} + {c16_8} + {c16_12} + {c16_16}) = {ans_a}")
print(f"({c12_4} + {c12_8} + {c12_12}) = {ans_c}")
print(f"(a) {ans_a}; (b) {ans_b}; (c) {ans_c}")