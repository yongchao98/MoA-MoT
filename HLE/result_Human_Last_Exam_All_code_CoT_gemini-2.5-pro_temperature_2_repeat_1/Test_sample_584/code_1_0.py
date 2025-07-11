import math

# Given parameters for the S-curve
r = 3
b = 9
false_positive_rate_limit = 0.05

print(f"We are given the S-curve P(s) = 1 - (1 - s^r)^b, with r={r} and b={b}.")
print(f"We need to find the similarity threshold 's' where the false positive rate is {false_positive_rate_limit}.")
print("-" * 60)

# Step 1: Set the S-curve formula equal to the false positive rate limit
print("Step 1: Set up the equation.")
print(f"1 - (1 - s^{r})^{b} = {false_positive_rate_limit}")

# Step 2: Rearrange to isolate the term containing 's'
one_minus_p = 1 - false_positive_rate_limit
print("\nStep 2: Isolate the parenthetical term.")
print(f"(1 - s^{r})^{b} = 1 - {false_positive_rate_limit}")
print(f"(1 - s^{r})^{b} = {one_minus_p}")

# Step 3: Take the b-th root of both sides
val1 = one_minus_p**(1 / b)
print("\nStep 3: Take the b-th root (1/b) of both sides.")
print(f"1 - s^{r} = ({one_minus_p})^(1/{b})")

# Step 4: Isolate the s^r term
val2 = 1 - val1
print(f"\nStep 4: Isolate the s^r term.")
print(f"s^{r} = 1 - ({one_minus_p})^(1/{b})")

# Step 5: Isolate s by taking the r-th root
s_threshold = val2**(1 / r)
print(f"\nStep 5: Isolate s by taking the r-th root (1/r).")
print(f"s = (1 - ({one_minus_p})^(1/{b}))^(1/{r})")

# Final calculation and result
print("\n--- Final Calculation ---")
print("The threshold 's' is calculated as:")
print(f"s = (1 - (1 - {false_positive_rate_limit})**(1/{b}))**(1/{r})")
print(f"s = {s_threshold:.3f}")

final_answer = f"{s_threshold:.3f}"
print(f"\nThe similarity threshold 's' should be {final_answer}.")
print(f"<<<{final_answer}>>>")