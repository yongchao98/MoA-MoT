import sys

def solve_euler_characteristic_mod_k():
    """
    This script computes the reduced Euler characteristic of a specially defined
    simplicial complex Delta_k modulo k, for a given prime k >= 3.
    """
    # We can get k from the command line, or use a default value.
    # The problem is stated for a general prime k, so we pick one for demonstration.
    k = 17 
    if len(sys.argv) > 1:
        try:
            k = int(sys.argv[1])
            if k < 3 or k % 2 == 0 and k != 2:
                print(f"Error: k must be an odd prime, but got {k}.")
                return
            # Simple primality test for odd numbers
            is_prime = all(k % i != 0 for i in range(3, int(k**0.5) + 1, 2))
            if not is_prime:
                 print(f"Error: k={k} is not a prime number.")
                 return
        except ValueError:
            print("Error: Please provide a valid integer for k.")
            return

    # Based on the reasoning using Smith's Fixed Point Theorem, the value of
    # the reduced Euler characteristic of Delta_k modulo k is (k-1)/2.

    # The number of non-empty fixed faces in the fixed-point subcomplex
    num_fixed_faces = (k - 1) // 2

    # The size of each fixed face is k
    face_size = k

    # The contribution of each fixed face to the Euler characteristic sum
    # is (-1)^(face_size - 1). Since k is an odd prime, k-1 is even.
    contribution_per_face = 1

    # The final result is the sum of contributions from all non-empty fixed faces.
    result = num_fixed_faces * contribution_per_face

    print(f"For the prime k = {k}:")
    print("The reduced Euler characteristic of Delta_k modulo k can be calculated by the formula:")
    print("result = (k - 1) / 2")
    print("Here is the step-by-step calculation:")
    print(f"result = ({k} - 1) / 2 = {k-1} / 2 = {result}")

solve_euler_characteristic_mod_k()