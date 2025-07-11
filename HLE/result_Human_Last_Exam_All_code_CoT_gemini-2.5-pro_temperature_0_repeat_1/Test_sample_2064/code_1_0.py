import math

def solve():
    """
    Calculates the final value based on the derived formula for l(a).

    The function l(a) is the derivative of the log-probability density function,
    l(a) = d/da [ln(p(U(a)))].
    Given the definitions, the log-PDF ln(p(Y)) can be expressed. For very large n,
    the dominant term in the derivative comes from the determinant of U(a).
    ln(p(U(a))) contains the term -ln(det(U(a))) = -ln(1 + n*a^2).
    The derivative of this term is l(a) = - (2*n*a) / (1 + n*a^2).
    For n = 1,000,000,000, this expression is extremely close to -2/a.
    l(a) approx - (2*n*a) / (n*a^2) = -2/a.
    The problem asks for floor(sum_{i=1 to 10} l(a_i)), where a_i are the first 10 prime numbers.
    We will calculate this sum using the approximation l(a_i) = -2/a_i.
    """
    
    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    # Calculate the sum of l(a_i)
    total_sum = 0
    
    # We print the equation step-by-step
    print("Calculating the sum S = floor(l(2) + l(3) + ... + l(29))")
    print("Using the approximation l(a) = -2/a for large n:")
    
    equation_parts = []
    for p in primes:
        term = -2.0 / p
        total_sum += term
        equation_parts.append(f"(-2/{p})")
        
    equation_str = "S = floor(" + " + ".join(equation_parts) + ")"
    print(equation_str)
    
    # Now calculate the numerical values and the sum
    sum_of_reciprocals = sum(1.0/p for p in primes)
    final_sum_val = -2 * sum_of_reciprocals
    
    print(f"S = floor(-2 * ({' + '.join([f'1/{p}' for p in primes])}))")
    print(f"S = floor(-2 * ({sum_of_reciprocals}))")
    print(f"S = floor({final_sum_val})")
    
    # Calculate the floor
    final_answer = math.floor(final_sum_val)
    
    print(f"The final answer is {final_answer}")
    
solve()