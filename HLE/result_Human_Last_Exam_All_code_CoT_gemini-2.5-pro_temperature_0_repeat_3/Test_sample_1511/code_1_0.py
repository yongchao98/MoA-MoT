import math

# Part (a): Based on the reasoning that a sum-maximal family can be F = all k-multisets
# and G = emptyset, and this F contains multisets with disjoint supports.
answer_a = "Yes"

# Part (b): Calculation for k=2, m=5.
# For m > k+1, the maximal sum is C(m+k-1, k).
m = 5
k = 2
m_plus_k_minus_1 = m + k - 1
result_b = math.comb(m_plus_k_minus_1, k)

# We print the calculation steps as requested.
print(f"Calculation for part (b):")
print(f"The maximal sum is given by the formula C(m+k-1, k).")
print(f"For m={m} and k={k}, this is C({m}+{k}-1, {k}) = C({m_plus_k_minus_1}, {k}).")
print(f"C({m_plus_k_minus_1}, {k}) = {result_b}")
answer_b = result_b

# Part (c): Based on the reasoning that the sum-maximal family F can be the empty set.
answer_c = "No"

# Print the final answer in the specified format.
print("\nFinal Answer:")
print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]")