import math

# Define the parameters for the S-curve
r = 3
b = 9

# Define the target false positive rate
p = 0.05

# To find the threshold, we need to solve the equation: 1 - (1 - s^r)^b = p
# Algebraically, s = (1 - (1 - p)^(1/b))^(1/r)

# Calculate the similarity threshold 's'
s_threshold = (1 - (1 - p)**(1/b))**(1/r)

# Print the final equation with each number, as requested
print(f"To find the threshold 's', we solve the S-curve equation with r={r}, b={b}, and a false positive rate of {p}:")
print(f"1 - (1 - s**{r})**{b} = {p}")

# Print the final calculated result rounded to three decimal points
print(f"\nThe similarity threshold 's' should be: {s_threshold:.3f}")