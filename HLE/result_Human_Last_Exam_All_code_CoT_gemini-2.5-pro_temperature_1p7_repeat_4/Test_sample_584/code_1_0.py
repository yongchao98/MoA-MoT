import math

# Given parameters for the S-curve
r = 3
b = 9
fp_rate_limit = 0.05

print(f"Given parameters: r = {r}, b = {b}")
print("The S-curve for the probability of a candidate pair is P(s) = 1 - (1 - s^r)^b")
print(f"We need to find the similarity threshold 's' where P(s) < {fp_rate_limit}.")
print("")
print("The inequality to solve is:")
print(f"1 - (1 - s**{r})**{b} < {fp_rate_limit}")
print("")
print("Solving for 's' gives us the expression:")
print(f"s < (1 - (1 - {fp_rate_limit})**(1/{b}))**(1/{r})")
print("")
print("Plugging in the numbers into the final equation step-by-step:")

# --- Calculation Steps ---
# Step 1: Calculate the base for the (1/b) exponent
val_step1 = 1 - fp_rate_limit
# Step 2: Calculate the result of the (1/b) exponent
val_step2 = val_step1**(1/b)
# Step 3: Calculate the base for the (1/r) exponent
val_step3 = 1 - val_step2
# Step 4: Calculate the final threshold value
threshold = val_step3**(1/r)

# --- Print Final Equation with Intermediate Values ---
print(f"s < (1 - ({val_step1})**(1/{b}))**(1/{r})")
print(f"s < (1 - {val_step2:.5f})**(1/{r})")
print(f"s < ({val_step3:.5f})**(1/{r})")
print(f"s < {threshold:.5f}")
print("")

# --- Final Answer ---
print(f"The threshold 's' must be less than {threshold:.3f} for the false positive rate to be below 0.05.")

final_answer = f"{threshold:.3f}"
print(f"<<<{final_answer}>>>")