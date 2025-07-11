import sys

def solve_euler_char_mod_k():
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k modulo k.

    The problem reduces to computing 4^{-1} mod k for a prime k >= 3.
    """
    print("This script computes hat_chi(Delta_k) mod k for a prime k >= 3.")
    
    try:
        k_str = input("Enter a prime k >= 3: ")
        k = int(k_str)
        
        if k < 3:
            print("Error: k must be greater than or equal to 3.")
            return

        # Basic primality test
        is_prime = True
        if k > 3:
            if k % 2 == 0:
                is_prime = False
            else:
                i = 3
                while i * i <= k:
                    if k % i == 0:
                        is_prime = False
                        break
                    i += 2
        
        if not is_prime:
            print(f"Error: {k} is not a prime number.")
            return
            
        # The result is 4^{-1} mod k.
        # This requires Python 3.8+ for the 3-argument pow with a negative exponent.
        # For older versions, one could use pow(4, k-2, k) based on Fermat's Little Theorem.
        if sys.version_info >= (3, 8):
            result = pow(4, -1, k)
        else:
            # Using Fermat's Little Theorem: a^(p-2) is the inverse of a mod p.
            result = pow(4, k - 2, k)

        print("\n--- Calculation ---")
        print(f"The reduced Euler characteristic of Delta_k modulo k is given by the formula:")
        print(f"hat_chi(Delta_{k}) mod {k} = -(k-1)/4 mod k = 1/4 mod k")
        
        print("\n--- Final Answer ---")
        print(f"For k = {k}, the result is:")
        print(f"hat_chi(Delta_{k}) mod {k} = {result}")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_euler_char_mod_k()
<<<1/4 mod k>>>