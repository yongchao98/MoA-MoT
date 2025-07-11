import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def get_order_of_root_of_unity(k, n):
    """Computes the order of the complex number exp(2*pi*i*k/n)."""
    if n == 0:
        return float('inf')
    common_divisor = gcd(k, n)
    return n // common_divisor

# --- Properties of the H3 Reflection Group ---
group_name = "H3"
group_order = 120
degrees = [2, 6, 10]
coxeter_number = max(degrees)

print(f"Analyzing the reflection group of type {group_name}")
print(f"Order of the group |H3|: {group_order}")
print(f"Coxeter number h: {coxeter_number}")
print(f"Degrees d_j: {degrees}")

# Exponents are d_j - 1
exponents = [d - 1 for d in degrees]
print(f"Exponents m_j: {exponents}\n")

print("Step 1: Identify candidate elements.")
print("An element must be 'regular' to have a 'regular eigenvector'.")
print("For an element to have an eigenvalue of order 10, its own order must be a multiple of 10.")
print("The orders of regular elements in H3 are divisors of h=10, so we only need to check regular elements of order 10.\n")

print("Step 2: Analyze the conjugacy classes of regular elements of order 10.")
print("In H3, there are two such classes: 10A (the Coxeter class) and 10B.\n")

# --- Analysis of Class 10A (Coxeter Class) ---
print("--- Analyzing Class 10A (Coxeter Class) ---")
print("The eigenvalues of a Coxeter element c are given by exp(2*pi*i * m_j / h).")

class_10A_eigenvalue_orders = [get_order_of_root_of_unity(m, coxeter_number) for m in exponents]
print(f"For exponents {exponents}, the orders of the eigenvalues are: {class_10A_eigenvalue_orders}")

has_order_10_eigenvalue_10A = 10 in class_10A_eigenvalue_orders
print(f"Does this class have an eigenvalue of order 10? {'Yes' if has_order_10_eigenvalue_10A else 'No'}\n")

# --- Analysis of Class 10B ---
print("--- Analyzing Class 10B ---")
print("Elements in this class are conjugate to c^3, where c is a Coxeter element.")
print("Their eigenvalues are the cubes of the eigenvalues of c.")

k = 3
class_10B_eigenvalue_orders = [get_order_of_root_of_unity(m * k, coxeter_number) for m in exponents]
print(f"The orders of the eigenvalues for class 10B are: {class_10B_eigenvalue_orders}")

has_order_10_eigenvalue_10B = 10 in class_10B_eigenvalue_orders
print(f"Does this class have an eigenvalue of order 10? {'Yes' if has_order_10_eigenvalue_10B else 'No'}\n")

print("Step 3: Conclude which elements satisfy the conditions.")
print("Only elements from the Coxeter class (10A) have eigenvalues of order 10.")
print("We must also confirm their corresponding eigenvectors are regular.")
print("An eigenvector for an eigenvalue exp(2*pi*i * m/h) is regular if gcd(m, h) = 1.")

is_condition_met = False
for m in exponents:
    order = get_order_of_root_of_unity(m, coxeter_number)
    is_regular = gcd(m, coxeter_number) == 1
    if order == 10 and is_regular:
        print(f"- For exponent m={m}, the eigenvalue has order {order} and its eigenvector is regular (gcd({m},{coxeter_number})=1).")
        is_condition_met = True

if is_condition_met:
    print("\nSince elements in class 10A meet the criteria, all elements in this class are the solutions.")
else:
    print("\nNo elements were found that satisfy the conditions.")

print("\nStep 4: Calculate the number of such elements.")
print("The number of elements is the size of the Coxeter conjugacy class (10A).")
print("Size = |H3| / h")

final_answer = group_order // coxeter_number

print(f"The final equation is: {group_order} / {coxeter_number} = {final_answer}")