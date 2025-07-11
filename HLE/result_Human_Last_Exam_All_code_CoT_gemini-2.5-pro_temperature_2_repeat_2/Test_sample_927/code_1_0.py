import decimal

def is_prime(num):
    """Checks if a number is prime."""
    if num < 2:
        return False
    for i in range(2, int(num**0.5) + 1):
        if num % i == 0:
            return False
    return True

def create_cs_parameter(s_set, prec=1000):
    """Creates the real parameter c_S for a given set S."""
    decimal.getcontext().prec = prec
    cs = decimal.Decimal(0)
    for k in s_set:
        term = decimal.Decimal(2)**(-3**k)
        cs += term
    return cs

def check_membership(n, cs):
    """
    Checks if n is in the set S by decoding the parameter c_S.
    This mimics the logic of the existential formula.
    Returns True if n is in S, False otherwise.
    """
    # In the formal proof, y = 2**(3**n) would be established
    # via an existential sub-formula for exponentiation.
    # Python's arbitrary-precision integers can compute this directly.
    y_val = decimal.Decimal(2)**(3**n) * cs

    # The floor function Z = floor(Y) is also existentially definable.
    z_val = int(y_val) # z_val is floor(y_val)

    # The parity check is existentially definable.
    is_odd = (z_val % 2 == 1)
    
    return is_odd

def main():
    """
    Main function to run the demonstration.
    """
    # 1. Define the set S. Let's use primes up to 20.
    s_set = {n for n in range(20) if is_prime(n)}
    print(f"The chosen set S is the set of prime numbers < 20: {sorted(list(s_set))}\n")

    # 2. Create the real parameter c_S for this set.
    # We need high precision because the terms 2**(-3**k) get very small.
    # Precision of 3**(max(S)+1) / log10(2) ~= 3**19 / 0.3 ~= 1000s of digits.
    # Let's set a sufficiently high precision for our example.
    prec = 3**(max(s_set) + 1)
    cs_param = create_cs_parameter(s_set, prec=prec)
    # print(f"The calculated parameter c_S is (first 100 digits): {str(cs_param)[:100]}...")

    # 3. Test membership for numbers from 0 to 19.
    print("Testing membership for n from 0 to 19:")
    print("-" * 40)
    print("n | Actual Member? | Decoded Member? | Result")
    print("-" * 40)
    for n_test in range(20):
        actual = (n_test in s_set)
        decoded = check_membership(n_test, cs_param)
        result = "OK" if actual == decoded else "FAIL"
        print(f"{n_test:<2}| {str(actual):<14} | {str(decoded):<15} | {result}")

if __name__ == '__main__':
    main()