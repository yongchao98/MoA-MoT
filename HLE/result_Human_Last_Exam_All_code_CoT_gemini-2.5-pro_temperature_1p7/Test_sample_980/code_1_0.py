def sum_proper_divisors(n):
    """
    Calculates the sum of the proper divisors of a number n.
    A proper divisor is a positive divisor of n other than n itself.
    """
    if n <= 1:
        return 0
    
    # Start with 1 as every number > 1 has 1 as a proper divisor
    total = 1
    divisors = [1]
    
    # Iterate from 2 up to the square root of n
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            # i is a divisor
            total += i
            divisors.append(i)
            # n/i is also a divisor
            pair = n // i
            if pair != i:
                total += pair
                divisors.append(pair)
    
    divisors.sort()
    return total, divisors

def solve_sequence():
    """
    Solves for the next number in the sequence based on the
    sum of proper divisors rule.
    """
    last_number = 124
    
    next_number, proper_divs = sum_proper_divisors(last_number)
    
    # Format the equation string
    equation_parts = [str(d) for d in proper_divs]
    equation_str = " + ".join(equation_parts)
    
    print(f"The rule is to find the sum of the proper divisors of the last number.")
    print(f"The last number is {last_number}.")
    print(f"The proper divisors of {last_number} are: {proper_divs}")
    print(f"The equation is: {equation_str} = {next_number}")
    print(f"Therefore, the next number in the sequence is {next_number}.")

solve_sequence()
<<<100>>>