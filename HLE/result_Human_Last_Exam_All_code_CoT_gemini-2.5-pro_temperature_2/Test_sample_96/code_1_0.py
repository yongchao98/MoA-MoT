# The number of torsion elements of order 10 in A/Z with minimal positive word length
# is equivalent to the number of elements of order 10 with minimal length in the Coxeter group W(E_8).

# Step 1: We determine the minimal length L.
# An element of order 10 can be formed by combining a 5-cycle and a 2-cycle.
# The minimal length of a 5-cycle is 4 (e.g., (1,2,3,4,5) = s_1*s_2*s_3*s_4 in S_5).
# The minimal length of a 2-cycle is 1 (a simple transposition, e.g., (6,7) = s_6 in S_8).
# If their supports are disjoint, their product has length 4 + 1 = 5.
# This occurs in the W(A_7) parabolic subgroup of W(E_8).
# So, the minimal length is 5.

# Step 2: We count the number of such elements.
# The elements are of the form w = w1 * w2, where:
# - w1 is a 5-cycle on {i, i+1, i+2, i+3, i+4} (or its inverse), which have length 4.
# - w2 is a 2-cycle on {j, j+1}, which has length 1.
# - The supports {i,...,i+4} and {j, j+1} must be disjoint.
# - The group is W(A_7) isomorphic to S_8, so the universe of indices is {1, ..., 8}.

total_elements = 0
case_results = []

print("We count the number of elements w in W(E_8) of minimal length and order 10.")
print("The minimal length of such an element is 5.")
print("These elements are products of a length-4 5-cycle and a length-1 2-cycle with disjoint supports on {1, ..., 8}.\n")

# Case 1: 5-cycle on {1,2,3,4,5}.
# There are 2 such cycles: (1,2,3,4,5) and its inverse (5,4,3,2,1).
# The disjoint 2-cycles on consecutive integers can be on {6,7} or {7,8}. So there are 2 choices.
count_case_1 = 2 * 2
total_elements += count_case_1
case_results.append(count_case_1)
print(f"Case 1: 5-cycle support is {{1,2,3,4,5}}. Choices for 2-cycle support: {{6,7}}, {{7,8}}.")
print(f"Number of elements = (2 choices for 5-cycle) * (2 choices for 2-cycle) = {count_case_1}")

# Case 2: 5-cycle on {2,3,4,5,6}.
# 2 cycle choices. Disjoint 2-cycles can only be on {7,8}. 1 choice.
count_case_2 = 2 * 1
total_elements += count_case_2
case_results.append(count_case_2)
print(f"Case 2: 5-cycle support is {{2,3,4,5,6}}. Choices for 2-cycle support: {{7,8}}.")
print(f"Number of elements = (2 choices for 5-cycle) * (1 choice for 2-cycle) = {count_case_2}")

# Case 3: 5-cycle on {3,4,5,6,7}.
# 2 cycle choices. Disjoint 2-cycles can only be on {1,2}. 1 choice.
count_case_3 = 2 * 1
total_elements += count_case_3
case_results.append(count_case_3)
print(f"Case 3: 5-cycle support is {{3,4,5,6,7}}. Choices for 2-cycle support: {{1,2}}.")
print(f"Number of elements = (2 choices for 5-cycle) * (1 choice for 2-cycle) = {count_case_3}")

# Case 4: 5-cycle on {4,5,6,7,8}.
# 2 cycle choices. Disjoint 2-cycles can be on {1,2} or {2,3}. 2 choices.
count_case_4 = 2 * 2
total_elements += count_case_4
case_results.append(count_case_4)
print(f"Case 4: 5-cycle support is {{4,5,6,7,8}}. Choices for 2-cycle support: {{1,2}}, {{2,3}}.")
print(f"Number of elements = (2 choices for 5-cycle) * (2 choices for 2-cycle) = {count_case_4}")

print(f"\nThe total number of such torsion elements is the sum of these cases.")
equation_str = " + ".join(map(str, case_results))
print(f"Total = {equation_str} = {total_elements}")
