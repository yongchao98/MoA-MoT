import sys

# Probabilities of edges BEING PRESENT
p_lower_h = 1.0  # Lower horizontal edge is always present
p_upper_h = 2/3  # Upper horizontal edge is present with p=2/3
p_vertical = 1/2   # Vertical edge is present with p=1/2

# In the limit c -> infinity, we calculate the expected time E[T] to advance
# one horizontal unit. The speed v = 1 / E[T].
# E[T] = pi_0 * t_0 + pi_1 * t_1
# where pi_y is the stationary probability of being at level y,
# and t_y is the expected time to advance starting from level y.

# As c->inf, t_0 -> 1, because the lower rail is always connected.
t_0 = 1

# The ratio pi_1 / pi_0 is found by balancing probability flow between levels.
# Upward flow rate from level 0 is proportional to p_vertical * exp(-c).
# Downward flow rate from level 1 is proportional to p_vertical * (p_upper_h*exp(-c) + (1-p_upper_h)*1).
# pi_0 * (p_vertical * exp(-c)) = pi_1 * (p_vertical * (1-p_upper_h))
# pi_1 / pi_0 = exp(-c) / (1-p_upper_h)
# However, this simplified model can be incorrect. A more careful calculation gives:
# E[P_up] = p_vertical * (1 / (e^c + e^{-c} + 1)) approx p_vertical * e^{-c}
# E[P_down] = p_vertical * E[1 / (I(h_left)e^{-c} + I(h_right)e^{c} + 1)]
# E[P_down] approx p_vertical * (p_upper_h * e^{-c} + (1 - p_upper_h) * 1)
# E[P_down] approx p_vertical * (1 - p_upper_h)
# pi_0 * p_vertical * e^{-c} = pi_1 * p_vertical * (1-p_upper_h)
# pi_1 / pi_0 = e^{-c} / (1-p_upper_h)
# The calculation in the thought block was pi_1/pi_0 = 3*exp(-c)
# Let's verify that calculation.
# E[P_{0->1}] approx (1/2) * e^{-c}
# E[P_{1->0}] approx (1/2) * ( (2/3)*e^{-c} + (1/3)*1 ) approx 1/6
# So pi_0 * (1/2) * e^{-c} = pi_1 * (1/6) => pi_1/pi_0 = 3 * e^{-c}
# This seems correct. Let's name the factor C1.
C1 = 3  # pi_1 / pi_0 approx C1 * exp(-c)

# t_1 is the expected time from level 1. It's dominated by traps.
# t_1 approx B * exp(c)
# A trap occurs when a walker at (n,1) moves to (n+1,1), but (n+1,1) is a dead end.
# This requires several edges to be present or absent.
# 1. Edge ((n,1),(n+1,1)) must exist. Prob = p_upper_h
# 2. Node (n+1,1) must be a dead end (no right, no down).
#    Edge ((n+1,1),(n+2,1)) must be missing. Prob = 1 - p_upper_h
#    Edge ((n+1,1),(n+1,0)) must be missing. Prob = 1 - p_vertical
#    Prob of dead end = (1 - p_upper_h) * (1 - p_vertical)
# 3. Escape must be possible from (n,1) for the loop to resolve.
#    Edge ((n,1),(n,0)) must exist. Prob = p_vertical
#
# Probability of this trapping configuration at (n,1) is:
# P(trap) = p_upper_h * p_vertical * (1 - p_upper_h) * (1 - p_vertical)
p_trap_config = p_upper_h * p_vertical * (1 - p_upper_h) * (1 - p_vertical)

# The time cost of this trap is approx 2 * exp(c).
time_cost_factor = 2

# So, B, the coefficient of exp(c) in t_1, is P(trap) * time_cost_factor.
B = p_trap_config * time_cost_factor

# The product pi_1 * t_1, in the limit c->inf, is:
# (C1 * exp(-c)) * (B * exp(c)) = C1 * B
limit_pi1_t1 = C1 * B

# And lim E[T] = lim (pi_0*t_0 + pi_1*t_1). As c->inf, pi_0 -> 1.
# lim E[T] = 1 * t_0 + limit_pi1_t1
limit_E_T = 1 * t_0 + limit_pi1_t1

# The speed v = 1 / E[T]
final_speed = 1 / limit_E_T

# We need to output the components of the final calculation
# Let's write the equation as numbers
# v = 1 / (1 + C1 * B) = 1 / (1 + 3 * (p_upper_h * p_vertical * (1-p_upper_h) * (1-p_vertical)) * 2)
# v = 1 / (1 + 3 * (2/3 * 1/2 * 1/3 * 1/2) * 2)
# v = 1 / (1 + 3 * (2/36) * 2)
# v = 1 / (1 + 3 * (1/18) * 2)
# v = 1 / (1 + (1/6) * 2)
# v = 1 / (1 + 1/3)
# v = 1 / (4/3)
# v = 3/4

print("The asymptotic speed is calculated as v = 1 / E[T], where E[T] is the expected time to advance one unit.")
print("For large c, E[T] approaches lim(pi_0*t_0 + pi_1*t_1).")
print("We found that t_0 = 1 and pi_0 approaches 1.")
print("The ratio pi_1/pi_0 approaches 3 * exp(-c).")
print("The time t_1 on the upper level is dominated by traps and is approx. (1/9) * exp(c).")
print("The product pi_1 * t_1 approaches (3 * exp(-c)) * ((1/9) * exp(c)) = 3/9 = 1/3.")
print("So, E[T] approaches 1 + 1/3 = 4/3.")
print("The final speed is v = 1 / (4/3) = 3/4.")
print("") # Blank line for spacing
print(f"Final equation: 1 / (1 + 3 * ({p_upper_h:.2f} * {p_vertical:.2f} * {1-p_upper_h:.2f} * {1-p_vertical:.2f}) * 2) = {final_speed}")
