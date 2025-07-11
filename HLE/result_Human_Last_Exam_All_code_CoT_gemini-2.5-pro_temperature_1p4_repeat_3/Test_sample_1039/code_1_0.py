import numpy as np
from fractions import Fraction

def calculate_variance_B3():
    """
    Calculates the variance of the Coxeter length statistic on the hyperoctahedral group B_3.

    The calculation is based on the length generating function for B_n:
    P(q) = Product_{i=1 to n} [2i]_q
    For B_3, this is P(q) = [2]_q * [4]_q * [6]_q.
    The coefficients of the expanded polynomial P(q) give the frequency of each length.
    """
    
    # Coefficients of the q-polynomials [2]_q, [4]_q, and [6]_q
    p2_coeffs = [1, 1]
    p4_coeffs = [1, 1, 1, 1]
    p6_coeffs = [1, 1, 1, 1, 1, 1]

    # Convolve the coefficient arrays to multiply the polynomials
    temp_coeffs = np.convolve(p2_coeffs, p4_coeffs)
    final_coeffs = np.convolve(temp_coeffs, p6_coeffs)

    # The final_coeffs array contains the number of elements for each length k at index k.
    # a_k = final_coeffs[k]
    
    num_elements = 0
    sum_of_lengths = 0
    sum_of_lengths_sq = 0

    print("Distribution of Coxeter lengths in B_3:")
    for length, count in enumerate(final_coeffs):
        print(f"Length {length}: {count} elements")
        num_elements += count
        sum_of_lengths += length * count
        sum_of_lengths_sq += length**2 * count
    
    print("\n--- Calculating Variance ---")
    
    # E[L] = (Σ l * a_l) / N
    mean = sum_of_lengths / num_elements
    
    # E[L^2] = (Σ l^2 * a_l) / N
    mean_sq = sum_of_lengths_sq / num_elements

    # Var(L) = E[L^2] - (E[L])^2
    variance = mean_sq - mean**2
    
    # For precise fractional representation
    mean_frac = Fraction(sum_of_lengths, num_elements)
    mean_sq_frac = Fraction(sum_of_lengths_sq, num_elements)
    variance_frac = mean_sq_frac - mean_frac**2

    print(f"\nTotal number of elements (N): {num_elements}")
    print(f"Sum of lengths (Σ l*a_l): {sum_of_lengths}")
    print(f"Sum of squared lengths (Σ l^2*a_l): {sum_of_lengths_sq}")

    print(f"\nMean length E[L] = {sum_of_lengths}/{num_elements} = {mean} (or {mean_frac})")
    print(f"Mean squared length E[L^2] = {sum_of_lengths_sq}/{num_elements} = {mean_sq:.4f} (or {mean_sq_frac})")
    
    print(f"\nVariance = E[L^2] - (E[L])^2")
    print(f"         = {mean_sq_frac} - ({mean_frac})^2")
    print(f"         = {mean_sq_frac} - {mean_frac**2}")
    print(f"         = {variance_frac}")
    print(f"         ≈ {variance:.10f}")


if __name__ == "__main__":
    calculate_variance_B3()