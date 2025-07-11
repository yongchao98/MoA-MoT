import numpy as np
from scipy.optimize import minimize_scalar

def solve():
    """
    This script verifies the lower bound of d(t,x) by analyzing a key polynomial inequality.
    """
    print("To determine a constant lower bound for d(t,x), we analyze the evolution of its minimum value.")
    print("This leads to a condition on the potential lower bound C, which must satisfy an inequality for all u in [0,1].")
    print("For the candidate bound C = -1, this inequality is K(u) >= 0, where:")
    print("K(u) = u^4 - u^3 - 5*u^2 + 3*u + 2\n")

    # Define the polynomial K(u)
    def K(u):
        return u**4 - u**3 - 5*u**2 + 3*u + 2

    # We verify that K(u) >= 0 for u in [0, 1] by finding its minimum value in this interval.
    # The analytical solution shows the minimum occurs at u=1.
    res = minimize_scalar(K, bounds=(0, 1), method='bounded')
    min_u = res.x
    min_val = res.fun

    print("We numerically find the minimum of K(u) on the interval [0, 1]:")
    print(f"The minimum of K(u) occurs at u = {min_u:.6f}")
    # We use a high precision format for the minimum value to show it's effectively zero.
    print(f"The minimum value is K({min_u:.6f}) = {min_val:.6e}\n")

    print("This numerical result supports the analytical finding that the minimum is 0.")
    print("Let's explicitly calculate the value of K(u) at the boundary u=1, which is the point of minimum.")
    
    u = 1
    # As requested, we output each number in the final equation calculation for K(1)
    print(f"Final equation check at u = {u}:")
    print(f"K({u}) = ({u})^4 - ({u})^3 - 5*({u})^2 + 3*({u}) + 2")
    
    val_at_1 = 1**4 - 1**3 - 5*(1)**2 + 3*(1) + 2
    print(f"     = 1 - 1 - 5 + 3 + 2")
    print(f"     = {val_at_1}\n")

    print("Since K(u) is non-negative for all u in [0, 1], C = -1 is a valid lower barrier for d(t,x).")
    print(f"The initial minimum value is d_min(0) = -0.5, which is greater than -1.")
    print("Therefore, by the minimum principle, d(t,x) must remain greater than or equal to -1 for all t >= 0.")
    
    lower_bound = -1.0
    print(f"\nThe constant lower bound of d(t,x) is {lower_bound}.")

solve()
<<< -1.0 >>>