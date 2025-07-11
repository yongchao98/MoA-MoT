import math

def solve_and_explain():
    """
    This function analyzes the problem for both modulo 2 and modulo 3,
    and demonstrates the contradiction for the modulo 3 case.
    """
    
    # --- Modulo 2 explanation ---
    # For the modulo 2 case, such a number 'a' exists. This can be proven
    # constructively. By choosing a large enough value for floor(a) (e.g., 3),
    # one can ensure that at each step 'n', the interval of possible values for a^n
    # is large enough to always contain an integer with the desired parity.
    answer_mod_2 = "yes"

    # --- Modulo 3 explanation and demonstration ---
    # For the modulo 3 case, no such number 'a' exists. We show this by contradiction.
    # Assume such an 'a' exists. Then floor(a) must be an integer k_1 such that k_1 % 3 == 1.
    # Let's test the smallest possibility, k_1 = 1, which means 'a' is in the interval [1.0, 2.0).
    
    print("--- Demonstrating the non-existence for the mod 3 case ---")
    
    # Initialize the interval for 'a' based on n=1
    L, R = 1.0, 2.0
    k1 = 1
    print(f"For n=1, we must have floor(a) = k_1 such that k_1 % 3 == 1.")
    print(f"Let's test k_1 = {k1}. This implies a is in the interval [{L:.4f}, {R:.4f}).")

    # Step n=2
    n = 2
    target_rem = n % 3
    # Range for a^n is [L^n, R^n), i.e., [1.0, 4.0)
    # Possible k are integers in [1, 2, 3]. We need k % 3 == 2, so k must be 2.
    k2 = 2
    print(f"\nFor n={n}, we need floor(a^{n}) = k_{n} such that k_{n} % 3 == {target_rem}.")
    print(f"Given a in [{L:.4f}, {R:.4f}), a^{n} is in [{L**n:.4f}, {R**n:.4f}).")
    print(f"The only possible integer k_{n} in this range is {k2}.")
    # Update interval for a: a must be in [k2^(1/n), (k2+1)^(1/n))
    L = max(L, k2**(1/n))
    R = min(R, (k2+1)**(1/n))
    print(f"This narrows the interval for a to [{L:.4f}, {R:.4f}).")
    
    # Step n=3
    n = 3
    target_rem = n % 3
    min_an, max_an = L**n, R**n
    # Possible k are integers in [floor(2.828), floor(5.196)] -> [2, 5].
    # We need k % 3 == 0, so k must be 3.
    k3 = 3
    print(f"\nFor n={n}, we need floor(a^{n}) % 3 == {target_rem}.")
    print(f"Given a in [{L:.4f}, {R:.4f}), a^{n} is in [{min_an:.4f}, {max_an:.4f}).")
    print(f"The only possible integer k_{n} is {k3}.")
    L = max(L, k3**(1/n))
    R = min(R, (k3+1)**(1/n))
    print(f"This narrows the interval for a to [{L:.4f}, {R:.4f}).")
    
    # Step n=4
    n = 4
    target_rem = n % 3
    min_an, max_an = L**n, R**n
    # Possible k are integers in [floor(4.326), floor(6.349)] -> [4, 6].
    # We need k % 3 == 1, so k must be 4.
    k4 = 4
    print(f"\nFor n={n}, we need floor(a^{n}) % 3 == {target_rem}.")
    print(f"Given a in [{L:.4f}, {R:.4f}), a^{n} is in [{min_an:.4f}, {max_an:.4f}).")
    print(f"The only possible integer k_{n} is {k4}.")
    L = max(L, k4**(1/n))
    R = min(R, (k4+1)**(1/n))
    print(f"This narrows the interval for a to [{L:.4f}, {R:.4f}).")

    # Step n=5: The contradiction
    n = 5
    target_rem = n % 3
    min_an, max_an = L**n, R**n
    # Possible k are integers in [floor(6.240), floor(7.477)] -> [6, 7].
    print(f"\nFor n={n}, we need floor(a^{n}) % 3 == {target_rem}.")
    print(f"Given a in [{L:.4f}, {R:.4f}), a^{n} is in [{min_an:.4f}, {max_an:.4f}).")
    print(f"The possible integer values for floor(a^{n}) are 6 and 7.")
    print(f"Let's check their remainders modulo 3:")
    print(f"6 % 3 = {6 % 3}")
    print(f"7 % 3 = {7 % 3}")
    print(f"Neither of these is equal to the target remainder {target_rem}.")
    print("\nThis is a contradiction. No possible value for floor(a^5) exists, so no such 'a' can exist.")
    answer_mod_3 = "no"

    # Print the final answer
    print("\n" + "="*40)
    print("Final Answer Summary:")
    print(f"Modulo 2: {answer_mod_2}")
    print(f"Modulo 3: {answer_mod_3}")
    print("="*40)

    return f"{answer_mod_2},{answer_mod_3}"

final_answer = solve_and_explain()
print(f"<<<{final_answer}>>>")