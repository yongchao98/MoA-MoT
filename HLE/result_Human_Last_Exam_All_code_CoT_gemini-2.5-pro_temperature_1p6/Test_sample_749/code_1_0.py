# Step 1: Define jump probabilities in the limiting (h=0, all-blue) environment.
p_L = 1/5  # Probability of jumping left
p_R = 4/5  # Probability of jumping right

print(f"In the h->0 limit (all-blue environment):")
print(f"Jump probability left, p_L = {p_L}")
print(f"Jump probability right, p_R = {p_R}")
print("-" * 20)

# Step 2: Calculate hitting probabilities for the sites adjacent to 0.
# The walk has a drift to the right, as p_R > p_L.
# For a particle starting at x > 0, the probability of ever hitting 0 is (p_L/p_R)^x.
rho_1 = (p_L / p_R)**1
# For a particle starting at x < 0, the drift is towards 0, so the probability of hitting 0 is 1.
rho_neg_1 = 1

print("Hitting probabilities for site 0:")
print(f"From site 1, rho_1 = ({p_L}/{p_R})^1 = {rho_1}")
print(f"From site -1, rho_-1 = {rho_neg_1}")
print("-" * 20)

# Step 3: Calculate the return probability to 0 for a particle starting at 0.
# The particle jumps to -1 with probability p_L or to 1 with probability p_R.
p_ret = p_L * rho_neg_1 + p_R * rho_1
print("The probability for a particle at 0 to return to 0 is:")
print(f"p_ret = p_L * rho_-1 + p_R * rho_1")
print(f"p_ret = {p_L} * {rho_neg_1} + {p_R} * {rho_1} = {p_ret}")
print("-" * 20)

# Step 4: Analyze the branching process of visitors to site 0.
# The mean number of offspring is m = h / (1 - p_ret).
# We take the limit as h -> 0.
# limit_m = lim_{h->0} h / (1 - p_ret)
# Since p_ret converges to a constant value 2/5, the limit is 0.
limit_m = 0

print("The mean offspring number 'm' for the branching process of visitors to site 0:")
print(f"m(h) = h / (1 - p_ret(h))")
print(f"As h -> 0, p_ret(h) -> {p_ret}")
print(f"So, lim_{h->0} m(h) = 0 / (1 - {p_ret}) = {limit_m}")
print("-" * 20)

# Step 5: Final conclusion.
# Since the mean offspring number is less than 1, the branching process is subcritical.
# The survival probability of a subcritical branching process is 0.
final_probability = 0
print(f"The process is subcritical, so the probability of infinite visitors is {final_probability}.")
