import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1_{P^n} tensored with O(2) on complex projective space P^n_C.

    The dimension 'n' of the projective space is not specified in the prompt,
    so we will use n=4 as an example. You can change this value to compute
    the dimension for any n >= 1.
    """
    n = 4

    # The dimension is given by the formula:
    # h^0(P^n, Omega^1(2)) = h^0(P^n, O(1)^{oplus(n+1)}) - h^0(P^n, O(2))

    # The first term, h^0(P^n, O(1)^{oplus(n+1)}), is (n+1)^2.
    # This corresponds to (n+1) copies of the space of linear homogeneous
    # polynomials in n+1 variables.
    term1 = (n + 1)**2

    # The second term, h^0(P^n, O(2)), is the dimension of the space of
    # quadratic homogeneous polynomials in n+1 variables, which is C(n+2, 2).
    # We use integer division // as the result must be an integer.
    term2 = (n + 2) * (n + 1) // 2
    
    # In Python 3.8+ you could also use math.comb:
    # term2 = math.comb(n + 2, 2)

    # The final result is the difference between these two terms.
    result = term1 - term2

    print(f"For n = {n}, the dimension of the space of global sections of Omega^1(2) is calculated as follows:")
    print(f"h^0(P^{n}, Omega^1(2)) = h^0(P^{n}, O(1)^{{oplus({n+1})}}) - h^0(P^{n}, O(2))")
    print("The final equation is:")
    print(f"{result} = {term1} - {term2}")

if __name__ == "__main__":
    calculate_dimension()
