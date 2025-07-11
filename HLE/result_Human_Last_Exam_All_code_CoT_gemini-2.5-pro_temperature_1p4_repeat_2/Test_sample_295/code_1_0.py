import math

def is_prime(n):
    """Helper function to check if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def compute_euler_char_mod_k():
    """
    Computes the reduced Euler characteristic of the complex Delta_k modulo k.
    """
    # Let k be a prime number >= 3. The user can change this value.
    k = 13

    if not (is_prime(k) and k >= 3):
        print(f"Error: k must be a prime number greater than or equal to 3.")
        print(f"The chosen k={k} is not valid.")
        return

    print(f"This program computes the reduced Euler characteristic of a simplicial complex Delta_k, modulo k, for k = {k}.")
    print("The method uses a fixed-point theorem from topological combinatorics.")
    print("-" * 20)

    # Step 1: Find the non-empty faces of the fixed subcomplex Delta_k^G.
    # These are the (k-1)/2 orbits of edges under the action of a k-cycle.
    num_faces = (k - 1) // 2
    print(f"The fixed subcomplex Delta_k^G has {num_faces} non-empty faces.")

    # Step 2: Each face corresponds to an orbit of size k.
    size_of_each_face = k
    print(f"Each of these faces consists of {size_of_each_face} edges.")

    # Step 3: Compute the reduced Euler characteristic of the fixed subcomplex, hat_chi(Delta_k^G).
    # This is given by the sum: sum (-1)^{|F|-1} for each face F.
    print("\nThe reduced Euler characteristic of the fixed subcomplex is the sum:")
    
    # We build and display the equation for hat_chi(Delta_k^G)
    term_str = f"(-1)^({size_of_each_face} - 1)"
    equation_list = [term_str] * num_faces
    
    print(f"hat_chi(Delta_k^G) = {' + '.join(equation_list)}")
    
    # Simplify the exponent
    term_val = (-1)**(size_of_each_face - 1)
    simplified_term_str = f"(-1)^({size_of_each_face - 1})"
    
    # Since k is an odd prime, k-1 is even, so (-1)^(k-1) = 1.
    print(f"Since k={k} is odd, the exponent {size_of_each_face - 1} is even, so {simplified_term_str} = {term_val}.")
    
    simplified_equation_list = [str(term_val)] * num_faces
    print(f"hat_chi(Delta_k^G) = {' + '.join(simplified_equation_list)}")

    # Calculate the final sum
    hat_chi_fixed_complex = num_faces * term_val
    print(f"hat_chi(Delta_k^G) = {hat_chi_fixed_complex}")
    print("-" * 20)

    # Step 4: Apply the theorem hat_chi(Delta_k) === hat_chi(Delta_k^G) (mod k).
    final_result_mod_k = hat_chi_fixed_complex % k
    
    print(f"By the theorem, hat_chi(Delta_k) is congruent to hat_chi(Delta_k^G) modulo {k}.")
    print(f"So, the result is {hat_chi_fixed_complex} mod {k}, which is {final_result_mod_k}.")
    print(f"The general formula for the result is (k-1)/2.")

compute_euler_char_mod_k()