import math

# Step 1: Define the jump probabilities and their ratios (rho)
p_left_red = 4/5
p_right_red = 1/5
rho_red = p_left_red / p_right_red

p_left_blue = 1/5
p_right_blue = 4/5
rho_blue = p_left_blue / p_right_blue

# The problem is for the limit h -> 0. We analyze the behavior for any small h > 0.
# As an example, let's take a representative small value for h.
h = 0.01 

# Step 2: Calculate the expected value of log(rho)
# E[log(rho)] = h * log(rho_red) + (1-h) * log(rho_blue)
# This simplifies to (2h - 1) * log(4)
expected_log_rho = (2 * h - 1) * math.log(rho_red)

# Step 3: Conclude based on the sign of the expectation
# If E[log(rho)] < 0, the random walk is transient to +infinity.
# If E[log(rho)] > 0, the random walk is transient to -infinity.
# If E[log(rho)] = 0, the random walk is recurrent.

# For any h in (0, 1/2), 2*h - 1 is negative, so E[log(rho)] is negative.
# This means the walk and the entire cloud of particles drift to +infinity.
# Therefore, site 0 can only be visited by a finite number of particles.
# The probability of infinitely many visits is 0 for any h in (0, 1/2).
# The limit as h -> 0 is therefore 0.

final_result = 0

print(f"The ratio of left-to-right jump probabilities at a red site (rho_red) is: {p_left_red} / {p_right_red} = {rho_red}")
print(f"The ratio of left-to-right jump probabilities at a blue site (rho_blue) is: {p_left_blue} / {p_right_blue} = {rho_blue}")
print(f"The condition for drift is based on E[log(rho)] = (2*h - 1) * log(4).")
print(f"For h={h}, E[log(rho)] = ({2*h - 1}) * log({rho_red}) = {expected_log_rho:.4f}, which is negative.")
print("Since this holds for all h in (0, 1/2), the process is always transient to +infinity.")
print("The probability of infinitely many particles visiting site 0 is 0.")
print(f"Therefore, the limit is: {final_result}")
