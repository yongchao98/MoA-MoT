from decimal import Decimal, getcontext

def solve_for_T():
    """
    Calculates the value of T based on the derived solvability condition.
    
    The numbers involved have extreme exponents, so we use Python's decimal
    module for arbitrary-precision arithmetic.
    """
    # Set the precision for the decimal calculations.
    # The magnitude is ~10^5010000, so we need enough precision to
    # represent the number accurately.
    # 5010000 digits for the integer part plus some for safety.
    getcontext().prec = 5010010

    # Define the given parameters as Decimal objects
    alpha = Decimal(10) ** 10000
    x0 = Decimal(10) ** -5000000
    y0 = Decimal(10) ** -5000000

    # The formula for T is:
    # T = alpha / [ (x0 / (1 - x0)) + (y0 / (1 - y0)) ]
    
    # Since x0 and y0 are equal, the denominator is 2 * (x0 / (1 - x0))
    sum_x = x0 / (Decimal(1) - x0)
    sum_y = y0 / (Decimal(1) - y0)
    
    denominator = sum_x + sum_y
    
    T = alpha / denominator
    
    print("Derived formula: T = alpha / [ (x0 / (1 - x0)) + (y0 / (1 - y0)) ]")
    print("\nValues used in the final equation:")
    # Printing such large numbers directly is not feasible. We represent them in scientific notation.
    print(f"alpha = 10^{10000}")
    print(f"x0 = 10^{-5000000}")
    print(f"y0 = 10^{-5000000}")
    
    # The result for T is 0.5 * 10^5010000 - 0.5 * 10^10000.
    # This number can be written as 4 followed by 5,000,000 nines, then a 5,
    # followed by 9,999 zeros. We print the scientific notation.
    print("\nCalculated value for T:")
    print("T is approximately 0.5 * 10^5010000")
    print(f"T = {T:e}")


solve_for_T()