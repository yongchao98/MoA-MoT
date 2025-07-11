import math

def get_p_part(n, p):
    """Calculates the p-part of an integer n (the highest power of p that divides n)."""
    if n == 0:
        return 1 # By convention
    
    p_power = 0
    while n > 0 and n % p == 0:
        p_power += 1
        n //= p
        
    return p**p_power

# Problem setup
# G = A_5 x C_2
order_A5 = math.factorial(5) // 2  # Order of Alternating group A_5
order_C2 = 2                       # Order of Cyclic group C_2
order_G = order_A5 * order_C2

# Characteristic of the field k
p = 2

# The dimension of the projective cover of the trivial module, dim(P_0),
# is the order of the Sylow p-subgroup of G.
dim_P0 = get_p_part(order_G, p)

# The principal block of kG has only one simple module, the trivial module L_0.
# The dimension of the trivial module is 1.
dim_L0 = 1

# The Cartan matrix is a 1x1 matrix (c_00), where c_00 = dim(P_0) / dim(L_0).
# The sum of all entries is simply c_00.
sum_of_entries = dim_P0 // dim_L0

print(f"The sum of all entries in the Cartan matrix is: {sum_of_entries}")
print("The final equation is derived as follows:")
print(f"{sum_of_entries} = {dim_P0} / {dim_L0}")