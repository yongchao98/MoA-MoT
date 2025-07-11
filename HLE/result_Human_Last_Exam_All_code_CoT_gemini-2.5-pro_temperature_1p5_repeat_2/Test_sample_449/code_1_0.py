import math

# Step 1: Define the starting point and calculate its distance from the origin.
x0_x = 3000
x0_y = 4000
norm_x0 = math.sqrt(x0_x**2 + x0_y**2)

# Step 2: Determine the constant C for the asymptotic expansion of the potential kernel a(x),
# using the normalization a(0)=0 and the exact value a(1,1) = 4/pi.
# a(x) is approx. (2/pi)log(||x||) + C
# C = a(1,1) - (2/pi)log(||(1,1)||) = 4/pi - (2/pi)log(sqrt(2)) = (4 - log(2))/pi
log_2 = math.log(2)
pi = math.pi
C = (4 - log_2) / pi

# Step 3: Calculate the value of the potential kernel at the starting point x0
# using the asymptotic formula.
log_norm_x0 = math.log(norm_x0)
a_x0 = (2 / pi) * log_norm_x0 + C

# Step 4: For the chosen normalization a(0)=0, we have a(1,0) = 1.
a_1 = 1

# Step 5: The probability of hitting the neighbors of the origin is a(1,0) / a(x0).
# The probability of never hitting is 1 - a(1,0) / a(x0).
prob_never_hit = 1 - (a_1 / a_x0)

# Print the results in a clear step-by-step fashion.
print(f"The starting point is ({x0_x}, {x0_y}), with distance from origin ||x0|| = {norm_x0:.0f}")
print(f"The potential kernel a(x) is normalized such that a(0)=0, which implies a(1,0) = {a_1}.")
print(f"The hitting probability is given by the formula P(hit A) = a(1,0) / a(||x0||).")
print(f"The desired probability is P(never hit A) = 1 - a(1,0) / a(||x0||).")
print(f"\nTo calculate a({norm_x0:.0f}), we use the asymptotic formula:")
print(f"a(x) ≈ (2/π) * log(||x||) + C")
print(f"where C = (4 - log(2))/π ≈ {C:.4f}")
print(f"a({norm_x0:.0f}) ≈ (2/π) * log({norm_x0:.0f}) + {C:.4f}")
print(f"a({norm_x0:.0f}) ≈ {a_x0:.4f}")
print(f"\nThe equation for the probability is:")
print(f"P(never hit A) = 1 - {a_1} / {a_x0:.4f}")
print(f"P(never hit A) ≈ {prob_never_hit:.4f}")
print(f"\nThe approximate answer with two significant digits is {prob_never_hit:.2f}.")
