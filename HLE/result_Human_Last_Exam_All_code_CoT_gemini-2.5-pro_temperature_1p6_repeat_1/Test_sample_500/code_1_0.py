import math

# Based on the step-by-step derivation, the long-term behavior of the number of customers
# in the system (X_t) depends on a parameter 'd' derived from the arrival rate (lambda)
# and the service time distribution.

# Given parameters:
lambda_rate = 3
# The service time tail probability for large u is P(S > u) = 1/(3*u) + m/(u*ln(u)).
# The problem states 'm' is a positive integer, so m >= 1.
# Let's consider the smallest possible value for m.
m = 1

# Step 1: Calculate the primary tail parameter 'c'.
# c = lim_{u->inf} u * P(S > u) = 1/3
c = 1/3

# Step 2: Check for criticality by computing lambda * c.
# If lambda * c = 1, the system is critical.
criticality_param = lambda_rate * c

# Step 3: For the critical case, analyze the second-order parameter 'd'.
# lambda * P(S>u) is of the form 1/u + d/(u*ln(u)).
# In our case, d = lambda * m.
d = lambda_rate * m

# Step 4: Apply the theorem for the liminf of X_t in the critical case.
# The theorem states that if d > 1/2, then liminf X_t = infinity.
# Let's perform this check.
threshold = 1/2

# The conclusion is that since m is a positive integer, d = 3*m will always be >= 3.
# Since 3 > 0.5, the condition d > 1/2 is always met.
# Therefore, the number of customers tends to infinity almost surely.
if d > threshold:
    result = "infinity"
else:
    # This case is not reachable under the problem's conditions.
    result = "Does not apply"

# The problem asks to calculate liminf X_t. Our derivation shows this is infinity.
# The 'final equation' is the inequality d > 1/2, where d = lambda * m.
# For m=1, the equation is 3*1 > 1/2, which is 3 > 0.5.
# The final answer is the result of this logical deduction.
print(result)