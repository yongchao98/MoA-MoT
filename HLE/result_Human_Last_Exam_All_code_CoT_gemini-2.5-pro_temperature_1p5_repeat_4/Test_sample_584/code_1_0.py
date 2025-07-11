import math

# --- Given Parameters ---
r = 3
b = 9
# This is the target probability (false positive rate) for a pair with similarity 's'
prob_s = 0.05

# --- The S-Curve Equation ---
# We need to find the similarity threshold 's' by solving the equation:
# prob_s = 1 - (1 - s^r)^b

# --- Algebraic Solution for s ---
# s = (1 - (1 - prob_s)^(1/b))^(1/r)

# --- Calculation ---
val1 = 1 - prob_s
val2 = val1**(1/b)
val3 = 1 - val2
s = val3**(1/r)

# --- Output the Process and Result ---
print("The threshold 's' is calculated by solving the S-curve equation for P(s) = 0.05:")
print(f"P(s) = 1 - (1 - s^r)^b")

# Here, we print the equation with the specific numbers filled in.
print(f"\nSolving the equation for s:")
print(f"{prob_s} = 1 - (1 - s^{r})^{b}")
print(f"(1 - s^{r})^{b} = 1 - {prob_s}")
print(f"(1 - s^{r})^{b} = {val1}")
print(f"1 - s^{r} = {val1}^(1/{b})")
print(f"1 - s^{r} = {val2:.5f}")
print(f"s^{r} = 1 - {val2:.5f}")
print(f"s^{r} = {val3:.5f}")
print(f"s = ({val3:.5f})^(1/{r})")
print(f"s = {s:.5f}")

# --- Final Answer ---
final_answer = round(s, 3)
print(f"\nThe threshold 's' where the false positive rate is < 0.05 should be: {final_answer}")
<<<0.179>>>