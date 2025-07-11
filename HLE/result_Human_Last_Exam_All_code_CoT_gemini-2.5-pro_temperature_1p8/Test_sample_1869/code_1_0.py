def solve():
    """
    Finds the smallest prime p > 3 in a list of choices such that
    Z[p-th root of 6] is not the ring of integers of Q(p-th root of 6).
    This condition is met if 6^(p-1) is congruent to 1 modulo p^2.
    """
    
    # Answer Choices
    primes = [17, 383, 1093, 66161, 534851]
    a = 6
    
    print(f"We are looking for the smallest prime p such that {a}^(p-1) == 1 (mod p^2).")
    print("-" * 30)

    for p in primes:
        p_squared = p * p
        
        # Calculate a^(p-1) mod p^2
        result = pow(a, p - 1, p_squared)
        
        print(f"Checking prime p = {p}:")
        print(f"The equation is {a}^({p}-1) mod {p}^2 = ?")
        print(f"The result is {result}")
        
        if result == 1:
            print(f"\nCondition met for p = {p}. This is the smallest such prime in the list.")
            final_equation = f"{a}^({p-1}) mod {p**2} = 1"
            print(f"The final equation is: {final_equation}")
            
            # This is the answer the user is looking for. Let's map it to the choices.
            if p == 17:
                print("The answer is A.")
            elif p == 383:
                print("The answer is B.")
            elif p == 1093:
                print("The answer is C.")
            elif p == 66161:
                print("The answer is D.")
            elif p == 534851:
                print("The answer is E.")

            return
        else:
            print("Condition not met.")
            print("-" * 30)

solve()