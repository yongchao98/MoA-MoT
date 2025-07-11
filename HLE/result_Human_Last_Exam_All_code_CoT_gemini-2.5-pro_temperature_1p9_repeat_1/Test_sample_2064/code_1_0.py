import math

def calculate_l_of_a(a):
    """
    Calculates the value of l(a) for the simplified case n=1.
    For n=1, the problem simplifies as follows:
    - M = [1], so Log(M) = [0], and the mean vector m is [0].
    - d = n(n+1)/2 = 1. Sigma is a 1x1 matrix.
    - Sigma_11 = (-1)^(1+1) * (min(1,1) - 1*1/(1+1)) = 1 * (1 - 1/2) = 1/2.
      So the distribution is Normal(0, 1/2).
    - U(a) = [1+a^2].
    - The random variable Y = Exp(X) where X ~ Normal(0, 1/2). Y is a log-normal variable.
    - The PDF of Y is p(y) = 1/(y * sqrt(2*pi*sigma^2)) * exp(-(ln(y))^2 / (2*sigma^2)).
      With sigma^2 = 1/2, p(y) = 1/(y*sqrt(pi)) * exp(-(ln(y))^2).
    - ln(p(y)) = -ln(y) - (ln(y))^2 + const.
    - We evaluate this at y(a) = 1+a^2.
    - l(a) = d/da [ln(p(U(a)))] = d/da [-ln(1+a^2) - (ln(1+a^2))^2].
    - Using the chain rule, d/da(ln(1+a^2)) = 2a/(1+a^2).
    - l(a) = -(2a/(1+a^2)) - 2*ln(1+a^2)*(2a/(1+a^2))
    - l(a) = -(2a/(1+a^2)) * (1 + 2*ln(1+a^2)).
    """
    term1 = 2 * a / (1 + a**2)
    term2 = 1 + 2 * math.log(1 + a**2)
    return -term1 * term2

def solve():
    """
    Solves the problem by assuming n=1, calculating l(a_i) for the first 10 primes,
    summing them, and taking the floor.
    """
    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    
    total_sum = 0
    
    # Header for the final equation output
    print("Calculating the sum: l(a_1) + l(a_2) + ... + l(a_10)")
    print("-" * 50)
    
    equation_parts = []
    for i, p in enumerate(primes):
        val = calculate_l_of_a(p)
        equation_parts.append(f"l({p}) = {val:.4f}")
        total_sum += val
    
    # Print each term in the final equation
    for part in equation_parts:
        print(part)
        
    print("-" * 50)
    # Print the sum of the equation
    print(f"Sum = {' + '.join([f'({p:.4f})' for p in [calculate_l_of_a(pr) for pr in primes]])}")
    print(f"Sum = {total_sum:.4f}")
    
    # Calculate the floor of the sum
    final_answer = math.floor(total_sum)
    
    print("-" * 50)
    print(f"Final calculation: floor(Sum) = floor({total_sum:.4f}) = {final_answer}")
    print("-" * 50)
    
    # As requested, output the final numerical answer.
    # print(f"The final answer is: {final_answer}")

solve()