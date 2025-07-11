import math

def calculate_phi_n():
    """
    This function solves the problem for a given n.
    The problem simplifies to calculating phi(n) = exp(2n - 4 + 2/n).
    The problem is stated for n >= 5. As no specific value of n is given,
    we will calculate the result for n=5 as a representative case.
    """
    n = 5
    
    # The final formula for the exponent is of the form a*n + b + c/n.
    # The coefficients are derived to be a=2, b=-4, c=2.
    # As requested, we output these numbers.
    a = 2
    b = -4
    c = 2
    print(f"The coefficients of the exponent formula a*n + b + c/n are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")

    # Calculate the exponent for n=5
    exponent = a * n + b + c / n
    
    # Calculate the final result phi(5)
    result = math.exp(exponent)
    
    print(f"\nFor n = {n}, the value of phi(n) is: {result}")

# Execute the function
calculate_phi_n()

# The final numerical answer for n=5
final_answer = math.exp(2 * 5 - 4 + 2 / 5)
# print(f"\n<<< {final_answer} >>>") # Suppressed for final output format