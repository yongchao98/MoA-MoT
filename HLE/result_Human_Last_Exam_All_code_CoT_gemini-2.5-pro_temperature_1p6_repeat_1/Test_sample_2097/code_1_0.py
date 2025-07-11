import math

def calculate_mz1_for_n(n):
    """
    Calculates the magnetization M_z(1) for a given number of spins n.

    The formula is derived by solving the Volterra integral equation:
    M_z(1, n) = [(-1)^n * 2^n / (n! * pi^n * n^n)] * S_n
    where S_n is the (n+1)-th finite difference of a factorial-like function.
    S_n = sum_{i=0}^{n+1} [ comb(n+1, i) * (-1)^(n+1-i) * (4n)! / (4n-i)! ]
    """
    # Calculate S_n using integer arithmetic to maintain precision
    s_n = 0
    # Python's math.factorial can handle large integers
    try:
        fact_4n = math.factorial(4 * n)
        for i in range(n + 2):
            sign = (-1)**(n + 1 - i)
            # This is (4n)! / (4n-i)! which is P(4n, i)
            # Use integer division // as the result is guaranteed to be an integer
            poly_term = fact_4n // math.factorial(4 * n - i)
            
            # math.comb also handles large numbers
            comb_term = math.comb(n + 1, i)
            s_n += comb_term * sign * poly_term
    except ValueError:
        # Handles cases where 4n-i < 0, although not expected in this loop
        return float('inf')

    # Calculate the prefactor components
    # Numerator coefficient of 1/pi^n
    num_coeff_val = (-1)**n * (2**n) * s_n
    
    # Denominator coefficient of 1/pi^n
    den_coeff_val = math.factorial(n) * (n**n)

    # Simplify the fraction by dividing by the greatest common divisor
    common_divisor = math.gcd(num_coeff_val, den_coeff_val)
    final_num = num_coeff_val // common_divisor
    final_den = den_coeff_val // common_divisor

    # Calculate the final float value
    mz_value = (final_num / final_den) / (math.pi**n)
    
    return mz_value, final_num, final_den

def find_minimum_magnetization():
    """
    Calculates M_z(1, n) for a range of n to find the minimum value.
    """
    min_mz = float('inf')
    n_min = -1
    min_details = {}
    
    # We test for n from 1 up to a reasonable limit (e.g., 10)
    # The function appears to have a single minimum in this range.
    for n_val in range(1, 11):
        mz_val, num, den = calculate_mz1_for_n(n_val)
        
        if mz_val < min_mz:
            min_mz = mz_val
            n_min = n_val
            min_details = {'n': n_min, 'num': num, 'den': den, 'val': min_mz}
            
    # Output the result as an equation
    print(f"The minimum magnetization occurs at n = {min_details['n']}.")
    
    num_part = min_details['num']
    den_part = min_details['den']
    n_part = min_details['n']
    val_part = min_details['val']
    
    print("The final equation for the minimum magnetization is:")
    print(f"M_z(1) = {num_part} / ({den_part} * pi^{n_part})")
    print("\nCalculated value:")
    print(f"M_z(1) = {val_part}")
    
    return val_part

if __name__ == '__main__':
    minimum_value = find_minimum_magnetization()
    # The final answer format as requested by the user prompt
    # The problem asks for the minimum magnetization value.
    # print(f"<<<{minimum_value}>>>")

find_minimum_magnetization()