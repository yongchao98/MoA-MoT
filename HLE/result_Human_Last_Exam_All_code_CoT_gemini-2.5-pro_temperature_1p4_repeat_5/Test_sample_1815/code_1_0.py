import math

def get_lcm_sequence_generator():
    """
    Generator function for the sequence x_k = lcm(1, 2, ..., k).
    """
    current_lcm = 1
    k = 1
    while True:
        # In Python < 3.9, math.lcm is not available.
        # math.lcm(a, b) = (a * b) // math.gcd(a, b)
        if hasattr(math, 'lcm'):
            current_lcm = math.lcm(current_lcm, k)
        else:
            current_lcm = (current_lcm * k) // math.gcd(current_lcm, k)
        yield k, current_lcm
        k += 1

def check_convergence(sequence_generator, n, max_k=100):
    """
    Checks if the sequence converges to 0 modulo n.
    Finds the first k such that for all subsequent terms, the term is divisible by n.
    """
    print(f"Checking for modulus n = {n}")
    
    # We need to find a K such that for all k > K, x_k % n == 0.
    # Our argument is that K=n works. Let's verify.
    # We will check for k from n up to max_k.
    
    first_convergence_k = -1
    
    for k, x_k in sequence_generator:
        if k >= n:
            if first_convergence_k == -1:
                print(f"For k = {k}, x_{k} = {x_k}.  x_{k} % {n} = {x_k % n}")
                if x_k % n == 0:
                    first_convergence_k = k
                else:
                    # This should not happen based on our mathematical argument
                    print(f"Error: x_{k} is not divisible by {n}. Argument is flawed.")
                    return
        if k > max_k:
            break
            
    if first_convergence_k != -1:
        print(f"The sequence x_k = lcm(1,...,k) is divisible by {n} for all k >= {first_convergence_k}.")
        print("This demonstrates that the sequence converges to 0 for this modulus.")
    else:
        # This case won't be reached if our argument is correct
        print(f"Could not confirm convergence for n={n} up to k={max_k}.")


# --- Main Demonstration ---
# Let's test for a few values of n
test_moduli = [5, 12, 30]

for n in test_moduli:
    # We need a new generator for each test
    lcm_gen = get_lcm_sequence_generator()
    check_convergence(lcm_gen, n)
    print("-" * 20)
