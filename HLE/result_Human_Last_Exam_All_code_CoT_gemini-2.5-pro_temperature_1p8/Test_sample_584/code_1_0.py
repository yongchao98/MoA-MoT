import math

# Parameters for the LSH S-curve
r = 3
b = 9
# The target false positive rate
p = 0.05

# The S-curve equation is P(s) = 1 - (1 - s^r)^b.
# We set this probability P(s) equal to our desired rate p and solve for the similarity s.
# The general solution for s is: s = (1 - (1 - p)**(1/b))**(1/r)

# Calculate the threshold s by substituting the given values
s_threshold = (1 - (1 - p)**(1/b))**(1/r)

# Round the result to three decimal points as requested
s_rounded = round(s_threshold, 3)

# Print the final equation with all the numerical values and the result.
# The first line shows the equation we need to solve.
# The second line shows the derived formula for 's' with all numbers substituted.
# The final line gives the calculated answer.
print(f"To find the threshold 's', we must solve the equation: 1 - (1 - s^{r})^{b} = {p}")
print(f"The expression for 's' is: s = (1 - (1 - {p})**(1/{b}))**(1/{r})")
print(f"The calculated threshold 's' rounded to three decimal points is: {s_rounded}")