import math

# --- Given Parameters ---
# r is the number of rows per band
r = 3
# b is the number of bands
b = 9
# fp_rate is the target false positive rate
fp_rate = 0.05

# --- Calculation ---
# The S-curve equation is P(s) = 1 - (1 - s^r)^b
# We need to solve for 's' when P(s) = fp_rate.
# The formula to solve for s is: s = (1 - (1 - fp_rate)^(1/b))^(1/r)
s_threshold = (1 - (1 - fp_rate)**(1/b))**(1/r)

# --- Output ---
print("The goal is to find the similarity threshold 's' where the S-curve probability is 0.05.")
print("The S-curve equation with r=3 and b=9 is:")
print(f"P(s) = 1 - (1 - s^{r})^{b}")
print("\nWe set the probability to the desired false positive rate and solve for 's':")
print(f"1 - (1 - s^{r})^{b} = {fp_rate}")

print("\nSolving for 's' gives the expression:")
print(f"s = (1 - (1 - {fp_rate})^(1/{b}))^(1/{r})")

print("\nFinal calculation:")
# Printing the equation with all numbers substituted as requested
print(f"s = (1 - ({1 - fp_rate})^(1/{b}))^(1/{r})")

# Printing the final answer rounded to three decimal places
print(f"\nThe calculated similarity threshold 's' is: {s_threshold:.3f}")
<<<0.179>>>