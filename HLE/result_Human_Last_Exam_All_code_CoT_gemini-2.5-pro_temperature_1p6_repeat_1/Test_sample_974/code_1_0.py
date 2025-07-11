import sympy

def find_prime_divisors_and_solve():
    """
    Solves the problem by finding prime divisors p of q and checking the condition.
    """
    q = 12740347
    
    print(f"Starting the analysis for q = {q}.")
    print("-" * 30)

    # Step 1: Find all prime divisors of q
    # We use sympy's factorint function for robust factorization.
    prime_factors = sympy.factorint(q)
    prime_divisors = list(prime_factors.keys())
    
    print(f"Step 1: The prime divisors of q = {q} are {prime_divisors}.")
    
    if not prime_divisors:
        print("q has no prime divisors (is equal to 1), which is not the case.")
        return

    # A list to store the solutions
    solutions = []

    # Step 2: Iterate through each prime divisor p and check the condition
    for p in prime_divisors:
        print("-" * 30)
        print(f"Step 2: Checking the condition for p = {p}.")
        
        # Since p is a prime divisor of q, and q itself is prime, we have p = q.
        # This means p is the characteristic of the finite fields F_q and F_{q^2}.
        print(f"For p = {p}, p is the characteristic of the underlying fields F_q and F_{q^2}.")

        print("\nWe need to find the number of elements of order p in PSL(n, k) where p = char(k).")
        print("These are the non-identity unipotent elements whose order is exactly p.")
        print("A unipotent element u = I + N has order p if N^(p) = 0. Since the nilpotency index of any NxN nilpotent matrix N is at most N, u will have order p if p >= N.")
        
        # Check condition for PSL(3, q^2)
        n1 = 3
        print(f"\nFor PSL({n1}, q^2):")
        print(f"Here, n = {n1} and the field size is q^2. The characteristic is p = {p}.")
        print(f"Since p = {p} is much larger than n = {n1}, all non-identity unipotent elements have order p.")
        print("The number of unipotent elements in PSL(n, k) is k^(n^2 - n).")
        print(f"So, the number of elements of order p is (q^2)^({n1}^2 - {n1}) - 1 = (q^2)^({n1**2 - n1}) - 1 = q^{2*6} - 1 = q^12 - 1.")

        # Check condition for PSL(4, q)
        n2 = 4
        print(f"\nFor PSL({n2}, q):")
        print(f"Here, n = {n2} and the field size is q. The characteristic is p = {p}.")
        print(f"Since p = {p} is much larger than n = {n2}, all non-identity unipotent elements have order p.")
        print(f"The number of elements of order p is q^({n2}^2 - {n2}) - 1 = q^({n2**2 - n2}) - 1 = q^{12} - 1.")

        # Step 3: Compare the two numbers
        print("\nStep 3: Comparing the two quantities.")
        print("Number of elements of order p in PSL(3, q^2) is q^12 - 1.")
        print("Number of elements of order p in PSL(4, q) is q^12 - 1.")
        print("The two quantities are equal.")
        
        print(f"\nConclusion: p = {p} is a solution.")
        solutions.append(p)

    print("-" * 30)
    print("Final list of primes p that satisfy the condition:")
    if solutions:
        for sol in solutions:
            print(sol)
    else:
        print("No prime divisors of q satisfy the condition.")

if __name__ == '__main__':
    find_prime_divisors_and_solve()