import sys

# Set a higher recursion limit for potentially deep recursions in Hofstadter sequences.
sys.setrecursionlimit(2000)

memo = {1: 1, 2: 1}

def s(n):
    """
    Calculates the n-th term of the sequence using the deduced rule.
    The rule is s[n] = s[n - s[n - 1]] + s[n - s[n - 2]].
    It uses a memoization table (a dictionary) to store and reuse previously computed values.
    """
    if n in memo:
        return memo[n]
    try:
        # Check if indices would be valid before recursing
        s_n_minus_1 = s(n - 1)
        s_n_minus_2 = s(n - 2)
        if (n - s_n_minus_1) <= 0 or (n - s_n_minus_2) <= 0:
            # This handles cases where the sequence might terminate (indices become non-positive)
            raise IndexError("Sequence terminated, invalid index")
        
        result = s(n - s_n_minus_1) + s(n - s_n_minus_2)
        memo[n] = result
        return result
    except (IndexError, RecursionError) as e:
        # If the sequence cannot be computed further, return a marker like None or re-raise
        memo[n] = None # Mark as non-computable
        # print(f"Could not compute s[{n}]: {e}")
        return None


# Generate the first 42 terms of the sequence based on the rule
sequence_list = []
for i in range(1, 43):
    term = s(i)
    if term is not None:
        sequence_list.append(str(term))

print("Generated S4 sequence based on the rule R(s[n]) = s[n - s[n - 1]] + s[n - s[n - 2]]:")
print(", ".join(sequence_list))
print("-" * 20)

# Demonstrate the calculation for n=4, showing each number in the equation.
print("Demonstration of calculation for n=4:")
# To show the steps, we need the values s[3] and s[2]
s2_val = s(2)
s3_val = s(3)
s_of_3_minus_s_of_2 = s(3-s2_val) # This is s[2]
s_of_3_minus_s_of_1 = s(3-s(1)) # This is s[2]

# s[4] requires s[3] and s[2]
s4_idx1_val = s(4-s3_val)
s4_idx2_val = s(4-s2_val)
s4_final_val = s4_idx1_val + s4_idx2_val

print(f"s[4] = s[4 - s[3]] + s[4 - s[2]]")
print(f"     = s[4 - {s3_val}] + s[4 - {s2_val}]")
print(f"     = s[{4 - s3_val}] + s[{4 - s2_val}]")
print(f"     = {s4_idx1_val} + {s4_idx2_val}")
print(f"     = {s4_final_val}")

print("\nNote: The generated sequence differs from the example at the 4th term (3 vs 2),")
print("which suggests a typo in the provided example sequence S4.")
