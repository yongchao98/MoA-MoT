import math

def calculate_sum_of_ratios():
    """
    Calculates the sum over all natural dimensions of the specified ratio.

    The problem defines a ratio R(d) for each dimension d. This function
    calculates R(d) for d=1, 2, 3, ... and sums them up. The series
    converges quickly, so we only need a finite number of terms for
    high precision.

    The final equation is S = R(1) + R(2) + R(3) + ...
    This code will print each number in this equation: the terms R(d) and the final sum S.
    """

    # Constant C_dist = E[sqrt(r_i^2 + r_j^2)] where r_i, r_j ~ U(-1, 1)
    # This is equal to (sqrt(2) + asinh(1)) / 3.
    c_dist = (math.sqrt(2) + math.asinh(1)) / 3

    total_sum = 0.0
    
    print("The final sum is the result of the infinite series S = R(1) + R(2) + R(3) + ...")
    print("The numbers in this equation are the terms R(d) and the final sum S.")
    print("-" * 20)

    # Loop over dimensions d. The terms decrease very rapidly, so d=20 is plenty.
    for d in range(1, 21):
        # Denominator of B(d) has a part (1 + (d-1)*C_dist)
        b_denom_part = 1 + (d - 1) * c_dist
        
        # Numerator of R(d) is (d+1)
        r_numerator = d + 1
        
        # Denominator of R(d) is (d! * 2^d * b_denom_part)
        try:
            r_denominator = math.factorial(d) * (2**d) * b_denom_part
        except OverflowError:
            # For large d, factorial and power will overflow, but the term is already ~0.
            r_d = 0.0
        else:
            r_d = r_numerator / r_denominator

        # Output each number (term) in the equation
        print(f"R({d}) = {r_d:.8f}")

        total_sum += r_d
        
        # Stop when terms become negligible to save computation
        if r_d < 1e-12:
            break
            
    print("-" * 20)
    # Output the final number in the equation (the sum S)
    print(f"The sum S converges to: {total_sum:.3f}")

if __name__ == '__main__':
    calculate_sum_of_ratios()