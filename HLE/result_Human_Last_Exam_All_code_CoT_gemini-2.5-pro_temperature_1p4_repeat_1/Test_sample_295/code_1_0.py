def solve_euler_characteristic_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of Delta_k modulo k for a prime k >= 3.

    Args:
        k: A prime number greater than or equal to 3.

    The problem asks for hat_chi(Delta_k) mod k. The simplicial complex Delta_k has
    the edges of the complete graph K_k as its ground set. A set of edges A is a face
    if the graph formed by these edges on the vertices of K_k has a maximum vertex
    degree of at most 2.

    The solution uses a result from equivariant algebraic topology. By considering the
    action of the cyclic group C_k on Delta_k, one can show that:
        hat_chi(Delta_k) = hat_chi(the subcomplex fixed by C_k) (mod k)

    The fixed subcomplex is analyzed and found to have a very simple structure, leading
    to the final formula:
        hat_chi(Delta_k) mod k = (k-1)/2

    This code implements this formula.
    """
    # As per the problem statement, k is a prime and k >= 3.
    # We assume the input k satisfies these conditions.

    numerator = k - 1
    denominator = 2

    # Since k is an odd prime, k-1 is even, so the division results in an integer.
    result = numerator // denominator

    print(f"For the prime k = {k}:")

    # The problem asks to output each number in the final equation.
    # The final equation is based on the derived formula: (k-1)/2
    print(f"The calculation is based on the derived formula (k-1)/2.")
    print(f"Substituting k = {k}, we get the final equation:")
    print(f"({k} - 1) / 2 = {result}")

    # The value (k-1)/2 is an integer between 0 and k-1, so it is the result modulo k.
    print(f"\nThe value of hat_chi(Delta_k) mod {k} is {result}.")


# Let's compute the result for a specific prime k.
# You can change this value to any other prime k >= 3.
k_value = 17
solve_euler_characteristic_mod_k(k_value)