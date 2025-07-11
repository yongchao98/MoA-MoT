import math

def get_primes_up_to(limit):
    """
    Generates a list of prime numbers up to a given limit using the
    Sieve of Eratosthenes.
    """
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[i]:
            for multiple in range(i * i, limit + 1, i):
                is_prime[multiple] = False
    
    primes_list = []
    for i in range(limit + 1):
        if is_prime[i]:
            primes_list.append(i)
    return primes_list

def solve_prime_sequence():
    """
    Calculates the 1st term of P^(11) by iteratively finding the
    n-th prime.
    """
    # The highest index for a prime we will need is x_10 = 648391.
    # p_n is approximated by n*ln(n). For n ~ 650,000, p_n is roughly 
    # 650000 * ln(650000) ~= 9.5 million. We use 11 million as a safe upper bound.
    sieve_limit = 11_000_000
    primes = get_primes_up_to(sieve_limit)
    
    # Let p(n) denote the n-th prime number.
    # Let x_k be the 1st term in P^(k).
    # The recurrence is x_{k+1} = p(x_k).

    print("Let p(n) denote the n-th prime number.")
    print("Let x_k denote the 1st term in the set P^(k).")
    print("-" * 30)
    
    # x_1 is the 1st prime number.
    k = 1
    # x_k stores the value of the term, which becomes the index for the next step.
    x_k_value = 2 
    print(f"The 1st term of P^({k}) is x_{k} = 2")

    for k in range(2, 12):
        prev_x_k_value = x_k_value
        
        # The new term's value is the prime at the index of the previous term's value.
        # Since the `primes` list is 0-indexed, p(n) is at primes[n-1].
        index = prev_x_k_value
        if index > len(primes):
            print(f"Error: The sieve limit was too low. Needed prime at index {index}.")
            return

        x_k_value = primes[index - 1]
        
        # Output the current step in the final equation format
        print(f"The 1st term of P^({k}) is x_{k} = p(x_{k-1}) = p({prev_x_k_value}) = {x_k_value}")

    print("-" * 30)
    print(f"The 1st term in P^(11) is {x_k_value}.")
    
# Execute the solver function
solve_prime_sequence()
<<<9836303>>>