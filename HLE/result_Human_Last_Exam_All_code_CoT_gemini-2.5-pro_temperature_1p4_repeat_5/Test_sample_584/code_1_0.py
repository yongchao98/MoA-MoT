import math

# --- Parameters ---
r = 3
b = 9
fp_rate = 0.05

# --- Explanation and Calculation ---
# The S-curve represents the probability P(s) that a pair of items with similarity 's'
# is flagged as a candidate. The formula is P(s) = 1 - (1 - s^r)^b.
# We are looking for the similarity threshold 's' where this probability
# is equal to the desired false positive rate.

print("Step-by-step calculation to find the similarity threshold 's':")
print("-" * 60)
print(f"Given parameters: r = {r}, b = {b}")
print(f"Target false positive rate (P_fp) = {fp_rate}")
print("\nWe solve the equation P(s) = P_fp for 's':")
print(f"1 - (1 - s^{r})^{b} = {fp_rate}")

print("\nStep 1: Isolate the term containing 's'")
print(f"(1 - s^{r})^{b} = 1 - {fp_rate}")
step1_val = 1 - fp_rate
print(f"(1 - s^{r})^{b} = {step1_val}")

print(f"\nStep 2: Take the {b}-th root of both sides")
print(f"1 - s^{r} = ({step1_val})^(1/{b})")
step2_val = step1_val**(1/b)
print(f"1 - s^{r} = {step2_val:.4f}")

print(f"\nStep 3: Isolate s^{r}")
print(f"s^{r} = 1 - {step2_val:.4f}")
step3_val = 1 - step2_val
print(f"s^{r} = {step3_val:.4f}")

print(f"\nStep 4: Take the {r}-rd root of both sides to find 's'")
print(f"s = ({step3_val:.4f})^(1/{r})")
final_s = step3_val**(1/r)
print(f"s = {final_s:.3f}")

print("-" * 60)
print(f"The threshold should be {final_s:.3f} to ensure the false positive rate is less than {fp_rate}.")
<<<0.179>>>