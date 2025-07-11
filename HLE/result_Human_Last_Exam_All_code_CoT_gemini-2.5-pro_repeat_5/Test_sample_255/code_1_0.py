import math

def euler_totient(n):
    """Computes Euler's totient function phi(n)."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

# Step 4: Calculate dimensions of H^1 groups
print("Step 4: Calculate dimensions of H^1 groups")

# For H^1(A,M) and H^1(B,M)
dim_H1_A = euler_totient(1)
print(f"The dimension of H^1(A, M) is phi(1) = {dim_H1_A}")

dim_H1_B = euler_totient(1)
print(f"The dimension of H^1(B, M) is phi(1) = {dim_H1_B}")

# For H^1(C,M)
divisors_of_8 = [d for d in range(1, 9) if 8 % d == 0]
dim_H1_C = 0
phi_values = []
for d in divisors_of_8:
    phi_d = euler_totient(d)
    phi_values.append(f"phi({d})={phi_d}")
    dim_H1_C += phi_d

print(f"The dimension of H^1(C, M) is the sum of phi(d) for d dividing 8:")
print(f"Divisors of 8 are: {divisors_of_8}")
print(f"Sum = {' + '.join(phi_values)} = {dim_H1_C}")
print(f"So, dim H^1(C, M) = {dim_H1_C}")

# Step 5: Determine the dimension of the image of alpha
print("\nStep 5: Determine the dimension of the image of alpha")
dim_im_alpha = dim_H1_A
print(f"The image of alpha is the subspace M_1, so its dimension is phi(1) = {dim_im_alpha}")

# Step 6: Final Calculation
print("\nStep 6: Final Calculation")
dim_H2_G_M = dim_H1_C - dim_im_alpha
print(f"The dimension of H^2(G, M) is dim H^1(C, M) - dim(im(alpha))")
print(f"Dimension = {dim_H1_C} - {dim_im_alpha} = {dim_H2_G_M}")
