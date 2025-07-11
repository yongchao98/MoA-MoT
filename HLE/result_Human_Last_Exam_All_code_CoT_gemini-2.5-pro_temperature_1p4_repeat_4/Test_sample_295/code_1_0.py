import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def compute_euler_characteristic_mod_k():
    """
    Computes the reduced Euler characteristic of the complex Delta_k, modulo k.

    The problem asks to compute chi_hat(Delta_k) mod k, where k >= 3 is a prime.
    The faces of Delta_k are edge sets of K_k forming graphs with max degree at most 2.
    Using a group action argument involving the cyclic group Z_k, it can be shown that:
    chi_hat(Delta_k) = (k-1)/2 (mod k).
    """
    try:
        k_str = input("Enter a prime number k >= 3: ")
        k = int(k_str)

        if k < 3 or not is_prime(k):
            print(f"Error: {k} is not a prime number greater than or equal to 3.")
            return

        # The derived formula for the result is (k-1)/2.
        # Since k is an odd prime, k-1 is even, so this is an integer division.
        numerator = k - 1
        denominator = 2
        result = numerator // denominator

        print(f"For k = {k}, the reduced Euler characteristic mod k is given by the formula (k-1)/2.")
        print("The final equation is:")
        print(f"({k} - 1) / {denominator} = {numerator} / {denominator} = {result}")
        print(f"So, chi_hat(Delta_{k}) mod {k} is {result}.")

    except ValueError:
        print("Invalid input. Please enter an integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    compute_euler_characteristic_mod_k()