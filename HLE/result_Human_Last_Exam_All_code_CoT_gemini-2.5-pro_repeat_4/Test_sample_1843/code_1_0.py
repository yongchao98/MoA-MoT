import math

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def print_step(title, content):
    """Helper function for printing formatted steps."""
    print(f"\n--- {title} ---")
    print(content)

# Step 1: Factorize N
N = 36036
factors_str = "2^2 * 3^2 * 7 * 11 * 13"
print_step("Step 1: Factorize the Conductor N",
           f"The conductor is N = {N}.\nIts prime factorization is N = {factors_str}.")

# Step 2: Theory of Primitive Characters
print_step("Step 2: Primitive Characters and Order",
           ("A Dirichlet character chi modulo N is primitive if it is a product of primitive characters\n"
            "chi = chi_4 * chi_9 * chi_7 * chi_11 * chi_13, where each chi_{p^k} is primitive for its modulus.\n"
            "The order of chi is lcm(ord(chi_4), ord(chi_9), ord(chi_7), ord(chi_11), ord(chi_13)).\n"
            "We need this lcm to be exactly 6."))

# Step 3: Count component character choices
print_step("Step 3: Count Choices for Each Component Character",
           ("We count the number of primitive characters for each prime-power modulus whose order divides k, for k=6, 3, and 2.\n"
            "Let n(p^k, d|k) be this number.\n"))

# Modulo 4
p4_ord2 = phi(2)
n_4_d6 = p4_ord2
n_4_d3 = 0
n_4_d2 = p4_ord2
print("Conductor 4 (2^2):")
print("  - The group of characters is cyclic of order phi(4)=2. There is 1 primitive character, which has order 2.")
print(f"  - Choices with order dividing 6: The order is 2, which divides 6. Count = {n_4_d6}.")
print(f"  - Choices with order dividing 3: The order is 2, which does not divide 3. Count = {n_4_d3}.")
print(f"  - Choices with order dividing 2: The order is 2, which divides 2. Count = {n_4_d2}.\n")

# Modulo 9
p9_ord3, p9_ord6 = phi(3), phi(6)
n_9_d6 = p9_ord3 + p9_ord6
n_9_d3 = p9_ord3
n_9_d2 = 0
print("Conductor 9 (3^2):")
print("  - The group of characters is cyclic of order phi(9)=6. Primitive characters have orders 3 and 6.")
print(f"  - Choices with order dividing 6: Orders 3 and 6 both divide 6. Count = phi(3)+phi(6) = {p9_ord3}+{p9_ord6} = {n_9_d6}.")
print(f"  - Choices with order dividing 3: Only order 3. Count = phi(3) = {n_9_d3}.")
print(f"  - Choices with order dividing 2: Neither 3 nor 6 divides 2. Count = {n_9_d2}.\n")

# Modulo 7
p7_ord2, p7_ord3, p7_ord6 = phi(2), phi(3), phi(6)
n_7_d6 = p7_ord2 + p7_ord3 + p7_ord6
n_7_d3 = p7_ord3
n_7_d2 = p7_ord2
print("Conductor 7:")
print("  - The group of characters is cyclic of order phi(7)=6. Primitive characters have orders 2, 3, and 6.")
print(f"  - Choices with order dividing 6: All primitive orders (2,3,6) divide 6. Count = phi(2)+phi(3)+phi(6) = {p7_ord2}+{p7_ord3}+{p7_ord6} = {n_7_d6}.")
print(f"  - Choices with order dividing 3: Only order 3. Count = phi(3) = {n_7_d3}.")
print(f"  - Choices with order dividing 2: Only order 2. Count = phi(2) = {n_7_d2}.\n")

# Modulo 11
p11_ord2 = phi(2)
n_11_d6 = p11_ord2
n_11_d3 = 0
n_11_d2 = p11_ord2
print("Conductor 11:")
print("  - The group of characters is cyclic of order phi(11)=10. Primitive orders are 2, 5, 10.")
print(f"  - Choices with order dividing 6: Only primitive order 2 divides 6. Count = phi(2) = {n_11_d6}.")
print(f"  - Choices with order dividing 3: No primitive order divides 3. Count = {n_11_d3}.")
print(f"  - Choices with order dividing 2: Only primitive order 2 divides 2. Count = phi(2) = {n_11_d2}.\n")

# Modulo 13
p13_ord2, p13_ord3, p13_ord6 = phi(2), phi(3), phi(6)
n_13_d6 = p13_ord2 + p13_ord3 + p13_ord6
n_13_d3 = p13_ord3
n_13_d2 = p13_ord2
print("Conductor 13:")
print("  - The group of characters is cyclic of order phi(13)=12. Primitive orders are 2, 3, 4, 6, 12.")
print(f"  - Choices with order dividing 6: Orders 2, 3, 6 divide 6. Count = phi(2)+phi(3)+phi(6) = {p13_ord2}+{p13_ord3}+{p13_ord6} = {n_13_d6}.")
print(f"  - Choices with order dividing 3: Only order 3. Count = phi(3) = {n_13_d3}.")
print(f"  - Choices with order dividing 2: Only order 2. Count = phi(2) = {n_13_d2}.\n")

# Step 4: Combine counts
print_step("Step 4: Calculate Total Counts using the Product Rule",
           ("The total number of primitive characters with conductor N and order dividing k, N(d|k),\n"
            "is the product of the component counts n(p^k, d|k)."))

N_d6 = n_4_d6 * n_9_d6 * n_7_d6 * n_11_d6 * n_13_d6
print(f"Number of characters with order dividing 6:")
print(f"  N(d|6) = {n_4_d6} * {n_9_d6} * {n_7_d6} * {n_11_d6} * {n_13_d6} = {N_d6}")

N_d3 = n_4_d3 * n_9_d3 * n_7_d3 * n_11_d3 * n_13_d3
print(f"Number of characters with order dividing 3:")
print(f"  N(d|3) = {n_4_d3} * {n_9_d3} * {n_7_d3} * {n_11_d3} * {n_13_d3} = {N_d3}")

N_d2 = n_4_d2 * n_9_d2 * n_7_d2 * n_11_d2 * n_13_d2
print(f"Number of characters with order dividing 2:")
print(f"  N(d|2) = {n_4_d2} * {n_9_d2} * {n_7_d2} * {n_11_d2} * {n_13_d2} = {N_d2}")

# Step 5: Final Calculation
print_step("Step 5: Final Calculation using Inclusion-Exclusion",
           ("The number of characters of order exactly 6, N(6), is given by:\n"
            "N(6) = N(d|6) - N(d|3) - N(d|2) + N(d|1)\n"
            "N(d|1) is 0, as primitive characters cannot have order 1."))

result = N_d6 - N_d3 - N_d2
print("\nThe final equation is:")
print(f"{N_d6} - {N_d3} - {N_d2} = {result}")
print("\nThus, the number of primitive Dirichlet characters of conductor 36036 and order 6 is 100.")
