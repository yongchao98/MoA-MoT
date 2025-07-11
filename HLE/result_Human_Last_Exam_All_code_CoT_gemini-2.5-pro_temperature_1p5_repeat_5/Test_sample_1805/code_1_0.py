# Set the integer n for which to calculate the expression nabla_q(T^n).
# You can change this value to any integer (e.g., 5, 1, 0, -2).
n = 5

# --- Main logic to derive and print the expression ---

print(f"For n = {n}:")

# Case 1: n = 0
# The q-derivative of a constant T^0=1 is 0.
if n == 0:
    # Equation: nabla_q(T^0) = 0
    print(f"nabla_q(T^{n}) = 0")

# Case 2: n = 1
# [1]_q = 1 and T^(1-1) = T^0 = 1. The result is 1.
elif n == 1:
    # Equation: nabla_q(T^1) = 1
    print(f"nabla_q(T^{n}) = 1")

# Case 3: n is a positive integer greater than 1
elif n > 1:
    # Build the q-integer polynomial: [n]_q = 1 + q + q^2 + ... + q^(n-1)
    q_terms = ["1"]
    for k in range(1, n):
        if k == 1:
            q_terms.append("q")
        else:
            q_terms.append(f"q^{k}")
    q_polynomial = f"({' + '.join(q_terms)})"

    # Build the T part: T^(n-1)
    power = n - 1
    if power == 1:
        t_part = "T"
    else:
        t_part = f"T^{power}"

    # Print the final equation, showing all numbers in the expression.
    # Example for n=5: nabla_q(T^5) = (1 + q + q^2 + q^3 + q^4) * T^4
    print(f"nabla_q(T^{n}) = {q_polynomial} * {t_part}")

# Case 4: n is a negative integer
else:  # n < 0
    m = -n  # m is a positive integer
    
    # The q-integer for negative n is [n]_q = -q^n * [-n]_q = -q^n * [m]_q
    # The sum form is -(q^-1 + q^-2 + ... + q^-m)
    q_terms = []
    for k in range(1, m + 1):
        q_terms.append(f"q^{-k}")
    q_laurent_polynomial = f"-({' + '.join(q_terms)})"

    # Build the T part: T^(n-1)
    power = n - 1
    t_part = f"T^{power}"

    # Print the final equation.
    # Example for n=-3: nabla_q(T^-3) = -(q^-1 + q^-2 + q^-3) * T^-4
    print(f"nabla_q(T^{n}) = {q_laurent_polynomial} * {t_part}")