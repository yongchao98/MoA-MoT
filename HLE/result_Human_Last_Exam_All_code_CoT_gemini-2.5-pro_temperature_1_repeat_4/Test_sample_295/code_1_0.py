import math

def is_prime(n):
    """Helper function to check if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def compute_reduced_euler_characteristic_mod_k(k):
    """
    Computes the reduced Euler characteristic of the complex Delta_k, modulo k.

    Let k >= 3 be a prime. We denote by K_k the complete graph on k vertices.
    Let Delta_k denote the abstract simplicial complex where the ground set is E(K_k)
    and a non-empty subset A of E(K_k) is a face if the graph (V(K_k), A)
    has all vertex degrees at most 2.

    The value to compute is chi_hat(Delta_k) = sum_{A in Delta_k, A non-empty} (-1)^{|A|-1} mod k.
    This quantity is equivalent to the standard Euler characteristic chi(Delta_k).

    We use a group theory argument. Let sigma be a k-cycle acting on the vertices of K_k.
    The Euler characteristic chi(Delta_k) is congruent to the Euler characteristic of the
    fixed-point subcomplex, chi(Delta_k^sigma), modulo k.

    A face A is fixed by sigma if and only if it is a union of edge orbits.
    The edge orbits are determined by the 'distance' between vertices on the cycle (1, 2, ..., k).
    There are (k-1)/2 such orbits, each of size k.
    A union of |J| such orbits forms a graph where every vertex has degree 2*|J|.
    For this to be a face, we need 2*|J| <= 2, which means |J| <= 1.

    So, the non-empty fixed faces are precisely the (k-1)/2 individual edge orbits.
    Each of these faces has size k, so its dimension is k-1.
    """
    print(f"Solving for k = {k}")
    
    # Validate input k
    if not isinstance(k, int) or k < 3 or not is_prime(k):
        print("Error: k must be a prime number greater than or equal to 3.")
        return

    # Step 1: Find the number of non-empty fixed faces.
    # These are the edge orbits. There are (k-1)/2 of them.
    num_fixed_faces = (k - 1) // 2
    print(f"Number of non-empty fixed faces = ({k} - 1) / 2 = {num_fixed_faces}")

    # Step 2: Find the dimension of these fixed faces.
    # Each fixed face is an orbit of size k, so its dimension is k-1.
    dim_fixed_faces = k - 1
    print(f"Dimension of each fixed face = {k} - 1 = {dim_fixed_faces}")

    # Step 3: Compute the Euler characteristic of the fixed-point subcomplex.
    # chi(Delta_k^sigma) = sum (-1)^dim * (number of faces with that dim)
    # Since all fixed faces have the same dimension, this is a single term.
    # Note: (-1)**(k-1) is 1 because k is an odd prime, so k-1 is even.
    chi_fixed_complex = ((-1)**dim_fixed_faces) * num_fixed_faces
    print(f"chi(Delta_k^sigma) = (-1)^({dim_fixed_faces}) * {num_fixed_faces} = {chi_fixed_complex}")

    # Step 4: The result is chi(Delta_k^sigma) mod k.
    result = chi_fixed_complex % k
    print(f"The value of chi_hat(Delta_k) mod {k} is equal to chi(Delta_k^sigma) mod {k}.")
    print(f"Final Result: {chi_fixed_complex} mod {k} = {result}")
    print("-" * 20)
    return result

if __name__ == '__main__':
    # Run the computation for a few example primes.
    compute_reduced_euler_characteristic_mod_k(3)
    compute_reduced_euler_characteristic_mod_k(5)
    compute_reduced_euler_characteristic_mod_k(17)
