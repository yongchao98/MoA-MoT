def power_adjoin_sqrt(base_u, base_v, d, exp):
    """
    Computes (base_u + base_v * sqrt(d))^exp.
    Returns a tuple (u, v) representing u + v*sqrt(d).
    """
    # (a + b*sqrt(d)) * (c + e*sqrt(d)) = (ac + bd*d) + (ae + bc)*sqrt(d)
    def multiply(u1, v1, u2, v2, d):
        new_u = u1 * u2 + v1 * v2 * d
        new_v = u1 * v2 + v1 * u2
        return new_u, new_v

    # Start with the multiplicative identity: 1 + 0*sqrt(d)
    res_u, res_v = 1, 0
    
    # Exponentiation by squaring
    while exp > 0:
        if exp % 2 == 1:
            res_u, res_v = multiply(res_u, res_v, base_u, base_v, d)
        base_u, base_v = multiply(base_u, base_v, base_u, base_v, d)
        exp //= 2
        
    return res_u, res_v

# We need to compute (3 + 1*sqrt(7))^20 = U_20 + V_20*sqrt(7)
u_base = 3
v_base = 1
d = 7
exponent = 20

u_20, v_20 = power_adjoin_sqrt(u_base, v_base, d, exponent)

# The sum of squares of coefficients is C_20 = U_20 + V_20
result = u_20 + v_20

print(f"Let (3 + sqrt(7))^20 = U_20 + V_20 * sqrt(7).")
print(f"U_20 = {u_20}")
print(f"V_20 = {v_20}")
print(f"The sum of the squares of the coefficients is U_20 + V_20.")
print(f"Result = {result}")