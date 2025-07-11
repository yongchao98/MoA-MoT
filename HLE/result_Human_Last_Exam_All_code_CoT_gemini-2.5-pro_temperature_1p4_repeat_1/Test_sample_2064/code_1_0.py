import math

def calculate_l_of_a(a):
    """
    Calculates the value of the function l(a) for a given a.
    This is based on the derived formula for the n=1 case.
    l(a) = -2a * (1 + 2*ln(1+a^2)) / (1+a^2)
    """
    a_squared = a * a
    log_val = math.log(1 + a_squared)
    numerator = -2 * a * (1 + 2 * log_val)
    denominator = 1 + a_squared
    return numerator / denominator

def solve():
    """
    Solves the problem by summing l(a_i) for the first 10 prime numbers and
    taking the floor of the result.
    """
    # First 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    total_sum = 0
    
    print("Calculating the sum of l(a_i) for the first 10 primes a_i:")
    # We print the intermediate values as requested by the prompt's hint
    # "Remember in the final code you still need to output each number in the final equation!"
    # The "final equation" is Sum = v1 + v2 + ... + v10
    equation_parts = []
    for p in primes:
        l_val = calculate_l_of_a(p)
        equation_parts.append(f"({l_val:.4f})")
        total_sum += l_val
    
    print("Sum = " + " + ".join(equation_parts))
    print(f"\nThe total sum is: {total_sum}")
    
    final_answer = math.floor(total_sum)
    print(f"The floor of the sum is: {final_answer}")
    
solve()