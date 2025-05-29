import numpy as np

# Given output
discount_factor = 0.7399828543832548
implied_rate = 0.0357072104407451

# Solve for T using the implied rate formula
T = -np.log(discount_factor) / implied_rate

# Solve for rate using the discount factor formula
rate = -np.log(discount_factor) / T

print({"rate": rate, "T": T})