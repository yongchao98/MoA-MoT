import decimal

def check_membership(n, c_A, precision=200):
    """
    Checks if a natural number n is in the set A encoded by c_A.

    The logic is that n is in A if and only if floor(2**(n+1) * c_A) is odd.
    """
    # Set precision for decimal arithmetic
    decimal.getcontext().prec = precision

    # Convert inputs to Decimal for high-precision arithmetic
    n = decimal.Decimal(n)
    
    # Calculate 2**(n+1)
    power_of_2 = decimal.Decimal(2) ** (n + 1)
    
    # Calculate the expression inside the floor function
    expression = power_of_2 * c_A
    
    # Calculate the floor of the expression
    floor_value = expression.to_integral_value(rounding=decimal.ROUND_FLOOR)
    
    # Check if the floor value is odd
    is_odd = (floor_value % 2 != 0)
    
    return is_odd

def main():
    """
    Demonstrates the encoding of an arbitrary set A into a real number c_A
    and the subsequent decoding to check for membership.
    """
    # Let's choose an arbitrary finite set of natural numbers for demonstration.
    # The argument holds for any subset of N, including infinite and non-computable sets.
    A = {0, 1, 4, 9, 16, 25} # The first few perfect squares
    print(f"The chosen set A is: {A}\n")
    
    # Set a high precision for the calculations
    # Needs to be greater than the largest element in A
    max_val_in_A = max(A) if A else 0
    precision = max_val_in_A + 50 
    decimal.getcontext().prec = precision

    # 1. Encode the set A into a single real number c_A
    c_A = decimal.Decimal(0)
    for num in A:
        term = decimal.Decimal(2) ** decimal.Decimal(-(num + 1))
        c_A += term

    print(f"The real parameter c_A (to {precision} digits precision) is approximately:")
    print(f"{c_A}\n")

    # 2. Decode c_A to check for membership
    print("Testing membership for numbers from 0 to 30:")
    print("=" * 40)
    print("n |  Is n in A? (decoded) | Actually in A?")
    print("-" * 40)
    
    for n_test in range(31):
        # The equation for checking membership: is floor(2**(n+1) * c_A) odd?
        # Here we simulate this check using our python function.
        is_member = check_membership(n_test, c_A, precision)
        
        # We need to output the numbers in the "final equation"
        power_of_2_val = decimal.Decimal(2) ** (decimal.Decimal(n_test) + 1)
        product_val = power_of_2_val * c_A
        floor_val = product_val.to_integral_value(rounding=decimal.ROUND_FLOOR)
        
        # Uncomment the line below for verbose output on each check
        # print(f"n={n_test}: floor(2**({n_test}+1) * c_A) = floor({product_val}) = {floor_val}")

        print(f"{n_test:<2}| {str(is_member):<23} | {str(n_test in A)}")
    print("="*40)
    
if __name__ == "__main__":
    main()
