import math

def encode_set(A):
    """
    Encodes a set of natural numbers A into a real number r_A.
    r_A = sum_{i in A} 2**(-i-1)
    """
    r_A = 0.0
    # For practical reasons with floating point precision, we only consider a finite subset
    # In theory, this works for any infinite subset as well.
    for i in A:
        r_A += 2**(-i-1)
    return r_A

def check_membership(i, r_A):
    """
    Checks if i is in the set encoded by r_A using the formula:
    b_i = floor(2**(i+1) * r_A) - 2 * floor(2**i * r_A)
    Returns True if b_i is 1, False otherwise.
    """
    val1 = 2**(i+1) * r_A
    k1 = math.floor(val1)
    
    val2 = 2**i * r_A
    k2 = math.floor(val2)
    
    # The bit b_i
    bit = k1 - 2 * k2
    
    # The final equation is k1 - 2*k2 = 1
    print(f"Checking for i={i}:")
    print(f"  k1 = floor(2**({i}+1) * r_A) = floor({val1:.4f}) = {k1}")
    print(f"  k2 = floor(2**{i} * r_A) = floor({val2:.4f}) = {k2}")
    print(f"  Final equation: {k1} - 2 * {k2} = {bit}")
    
    return bit == 1

def main():
    # Let's define an arbitrary set of natural numbers
    # For this example, let's take the set of prime numbers less than 20
    A = {2, 3, 5, 7, 11, 13, 17, 19}
    print(f"The original set A is: {A}\n")

    # Encode this set into a real number parameter r_A
    r_A = encode_set(A)
    print(f"The real number parameter r_A encoding A is approximately: {r_A}\n")

    # Now, let's check for membership of a few numbers
    # using only r_A and the decoding formula
    test_numbers = [2, 4, 5, 20]
    for i in test_numbers:
        is_member = check_membership(i, r_A)
        print(f"--> Is {i} in A? {is_member}\n")

if __name__ == "__main__":
    main()
