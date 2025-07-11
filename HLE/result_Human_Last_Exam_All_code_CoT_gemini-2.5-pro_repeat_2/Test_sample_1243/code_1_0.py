import math

# Part (a)
# We count the number of non-zero binary vectors v in (Z/2Z)^16 whose weight w(v) is a multiple of 4.
# The possible weights are k = 4, 8, 12, 16.
# The number is C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16).
c_16_4 = math.comb(16, 4)
c_16_8 = math.comb(16, 8)
c_16_12 = math.comb(16, 12)
c_16_16 = math.comb(16, 16)
a_val = c_16_4 + c_16_8 + c_16_12 + c_16_16

# Part (b)
# This is a theoretical question. As explained in the plan, such a neighbor cannot be even.
b_val = "no"

# Part (c)
# We count the number of binary vectors v in (Z/2Z)^12 with weight 4.
# The number is C(12, 4).
c_val = math.comb(12, 4)

# Print the final answer in the required format.
print(f"(a) The number of distinct 2-isotropic vectors is given by the sum of counts for vectors of weight 4, 8, 12, and 16.")
print(f"C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16) = {c_16_4} + {c_16_8} + {c_16_12} + {c_16_16} = {a_val}")

print(f"(b) As reasoned in the plan, it is not possible for the neighbor to be even. The answer is {b_val}.")

print(f"(c) The number of unimodular 2-neighbors corresponds to the number of characteristic vectors of weight 4 in dimension 12.")
print(f"C(12, 4) = {c_val}")

# Final answer in the specified format
final_answer = f"(a) {a_val}; (b) {b_val}; (c) {c_val}"
print(f"\n<<<{final_answer}>>>")