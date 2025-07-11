from fractions import Fraction

# Step 1: Explain the model and the method
print("This problem involves three particles in a continuous-time random walk with an exclusion interaction.")
print("The system's state is described by the gaps between adjacent particles: Y1 = X2 - X1 and Y2 = X3 - X2.")
print("The stationary distribution for these gaps is a product of two geometric distributions.")
print("P(Y1=k) is proportional to alpha^(k-1) and P(Y2=k) is proportional to beta^(k-1).")
print("We find alpha and beta by solving the steady-state zero-flux balance equations.\n")

# Step 2: Set up and solve the linear equations for alpha and beta
print("The balance equations for the parameters alpha and beta form a system of two linear equations:")
print("1) 2*alpha - beta = 1/3")
print("2) -alpha + 2*beta = 1\n")

# We solve this system symbolically.
# From (2), we get alpha = 2*beta - 1.
# Substitute into (1): 2*(2*beta - 1) - beta = 1/3
# 4*beta - 2 - beta = 1/3  =>  3*beta = 7/3  =>  beta = 7/9
# Then, alpha = 2*(7/9) - 1 = 14/9 - 9/9 = 5/9.
alpha = Fraction(5, 9)
beta = Fraction(7, 9)

print("Solving the system of equations yields:")
print(f"alpha = {alpha}")
print(f"beta = {beta}\n")

# Step 3: Calculate the average gaps
# For a geometric distribution P(k) = p*(1-p)^(k-1) for k>=1, the mean is 1/p.
# Here, the parameter p is (1-alpha) for Y1 and (1-beta) for Y2.
p_Y1 = 1 - alpha
E_Y1 = 1 / p_Y1

p_Y2 = 1 - beta
E_Y2 = 1 / p_Y2

print("The average gaps are the means of their respective geometric distributions.")
print(f"The probability of the first gap being 1 is p1 = 1 - alpha = 1 - {alpha} = {p_Y1}.")
print(f"The average first gap is E[Y1] = 1 / p1 = 1 / {p_Y1} = {E_Y1}\n")

print(f"The probability of the second gap being 1 is p2 = 1 - beta = 1 - {beta} = {p_Y2}.")
print(f"The average second gap is E[Y2] = 1 / p2 = 1 / {p_Y2} = {E_Y2}\n")

# Step 4: Calculate the total average distance
total_distance = E_Y1 + E_Y2
print("The average distance between the leftmost and the rightmost particles is E[Y1] + E[Y2].")
print(f"Average Distance = {E_Y1} + {E_Y2} = {total_distance}\n")

# Step 5: Calculate the asymptotic speed of the leftmost particle
# The speed is v1 = (rate of jumping right) * (+1) + (rate of jumping left) * (-1)
# The jump to the right is suppressed if Y1 = 1.
# v1 = (1 * P(Y1 > 1)) - (1/3)
# P(Y1 > 1) = 1 - P(Y1=1) = 1 - p_Y1 = alpha.
rate_right = alpha
rate_left = Fraction(1, 3)
speed = rate_right - rate_left
print("The asymptotic speed of the leftmost particle is calculated from its effective jump rates.")
print(f"v1 = (rate right) - (rate left) = (1 * P(Y1 > 1)) - 1/3")
print(f"The probability of a right jump (Y1 > 1) is P(Y1 > 1) = alpha = {alpha}")
print(f"v1 = {alpha} - {rate_left} = {speed}\n")

# Step 6: Final Answer
print("The final calculated results are:")
print(f"Average distance = {total_distance}")
print(f"Asymptotic speed = {speed}")