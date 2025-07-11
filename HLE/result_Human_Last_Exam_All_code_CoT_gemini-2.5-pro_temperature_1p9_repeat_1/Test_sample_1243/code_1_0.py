import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    """
    if k < 0 or k > n:
        return 0
    # math.comb is efficient and handles large numbers
    return math.comb(n, k)

# Part (a)
# For n=16, d=2, we need to count non-zero vectors v in {0,1}^16
# where weight(v) is a multiple of 4.
# Possible weights are 4, 8, 12, 16.
n_a = 16
weights_a = [4, 8, 12, 16]
counts_a = [combinations(n_a, k) for k in weights_a]
result_a = sum(counts_a)

# Part (b)
# The reasoning is theoretical and the conclusion is definitive.
result_b = "no"

# Part (c)
# For n=12, d=2, we need to count non-zero vectors v in {0,1}^12
# where weight(v) is a multiple of 4.
# Possible weights are 4, 8, 12.
n_c = 12
weights_c = [4, 8, 12]
counts_c = [combinations(n_c, k) for k in weights_c]
result_c = sum(counts_c)

# --- Output the results ---
# Output the step-by-step calculation for each part as requested.
print("Calculation Details:")

# Part (a) details
eq_str_a = " + ".join([f"C(16, {k})" for k in weights_a])
val_str_a = " + ".join(map(str, counts_a))
print(f"(a) The number is {eq_str_a} = {val_str_a} = {result_a}")

# Part (b) details
print(f"(b) A theoretical proof shows this is impossible. The answer is '{result_b}'.")

# Part (c) details
eq_str_c = " + ".join([f"C(12, {k})" for k in weights_c])
val_str_c = " + ".join(map(str, counts_c))
print(f"(c) The number is {eq_str_c} = {val_str_c} = {result_c}")

# Output the final answer in the specified format.
print("\n---")
print(f"(a) [{result_a}]; (b) [{result_b}]; (c) [{result_c}]")