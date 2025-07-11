import numpy as np
from scipy.special import gamma

def calculate_values_for_case2():
    """
    Calculates the LHS and a lower bound for the RHS of the inequality for function 2.
    """
    # K1 = Integral[0,1] d(xi)/sqrt(xi(1-xi^2)) = Gamma(1/4)^2 / (2*sqrt(2*pi))
    K1 = (gamma(1/4)**2) / (2 * np.sqrt(2 * np.pi))

    # LHS = Area / pi = (K1^2 / 2) / pi
    lhs = (K1**2) / (2 * np.pi)

    # a0 = f(0), |a0| = K1/4
    abs_a0 = K1 / 4
    
    # a1 = f'(0), |a1| = sqrt(2)
    abs_a1 = np.sqrt(2)
    
    # RHS is sum(|a_n|), which has a lower bound of |a0| + |a1|
    rhs_lower_bound = abs_a0 + abs_a1
    
    print("Analysis for Function 2:")
    print(f"The constant K1 is approximately: {K1:.4f}")
    print(f"LHS = Area/pi = K1^2/(2*pi) is approximately: {lhs:.4f}")
    print(f"Lower bound for RHS >= |a0| + |a1| is approximately: {rhs_lower_bound:.4f}")
    
    print("\nFinal inequality check for Function 2:")
    print(f"Is {lhs:.4f} <= {rhs_lower_bound:.4f}? {lhs <= rhs_lower_bound}")
    print("The full sum on the RHS is even larger, so the inequality holds for function 2.")
    
    # The final equation is LHS <= RHS. Printing the components.
    print(f"\nThe equation is sum(n*|a_n|^2) <= sum(|a_n|).")
    print("We compared the following values to show this holds:")
    print(f"{lhs} <= {rhs_lower_bound} + ...")

calculate_values_for_case2()