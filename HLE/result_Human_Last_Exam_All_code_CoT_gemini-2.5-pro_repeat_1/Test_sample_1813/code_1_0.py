import math

def get_rational_cf(p, q):
    """Computes the continued fraction of a rational number p/q."""
    coeffs = []
    # Ensure p < q for the form [0; a1, a2, ...]
    if p >= q:
        raise ValueError("This function assumes p < q for simplicity.")
    
    # First coefficient is always 0
    coeffs.append(p // q)
    p = p % q
    
    while q != 0:
        a = p // q
        # This is an intermediate step, the actual coefficient is from q/p
        p, q = q, p
        a = p // q
        coeffs.append(a)
        p = p % q
        
    return coeffs

def get_quadratic_equation_from_period(period):
    """
    Computes the quadratic equation Ax^2 + Bx + C = 0 for the value x of a
    purely periodic continued fraction [0; period].
    The value of the fraction x is defined as x = [0; period].
    Let theta = [period; period]. Then x = 1/theta.
    """
    # We use the matrix method to find the equation for theta.
    # For a continued fraction [c1, c2, ...], the matrix is M1*M2*...
    # where Mi = [[ci, 1], [1, 0]].
    # The final matrix is [[p_k, p_{k-1}], [q_k, q_{k-1}]].
    total_m = [[1, 0], [0, 1]]  # Start with identity matrix

    for c_i in period:
        m_i = [[c_i, 1], [1, 0]]
        # Multiply total_m by m_i: total_m = total_m * m_i
        a, b = total_m[0]
        c, d = total_m[1]
        next_a = a * m_i[0][0] + b * m_i[1][0]
        next_b = a * m_i[0][1] + b * m_i[1][1]
        next_c = c * m_i[0][0] + d * m_i[1][0]
        next_d = c * m_i[0][1] + d * m_i[1][1]
        total_m = [[next_a, next_b], [next_c, next_d]]

    p_k, p_km1 = total_m[0]
    q_k, q_km1 = total_m[1]

    # The value theta = [period; period] satisfies the equation:
    # q_k * theta^2 + (q_km1 - p_k) * theta - p_km1 = 0
    A_theta = q_k
    B_theta = q_km1 - p_k
    C_theta = -p_km1

    # The value of our continued fraction is x = [0; period] = 1/theta.
    # Substituting theta = 1/x into the equation gives:
    # A_theta/x^2 + B_theta/x + C_theta = 0
    # Multiplying by x^2 gives: A_theta + B_theta*x + C_theta*x^2 = 0
    # Rearranging gives the standard form: C_theta*x^2 + B_theta*x + A_theta = 0
    A_x = C_theta
    B_x = B_theta
    C_x = A_theta

    # Simplify the equation by dividing by the GCD of coefficients.
    common_divisor = math.gcd(math.gcd(A_x, B_x), C_x)
    A_x //= common_divisor
    B_x //= common_divisor
    C_x //= common_divisor

    # Ensure the leading coefficient is positive for a standard representation.
    if A_x < 0:
        A_x, B_x, C_x = -A_x, -B_x, -C_x
        
    return A_x, B_x, C_x

# --- Main Execution ---

# 1. Define the rational number
p, q = 4, 7

# 2. Find its continued fraction: [0; 1, 1, 3]
cf_rational_coeffs = get_rational_cf(p, q)
rational_part = cf_rational_coeffs[1:]
print(f"The continued fraction of {p}/{q} is [{cf_rational_coeffs[0]}; {', '.join(map(str, rational_part))}].")

# 3. Construct the periodic part of the associated quadratic irrationality
# The period is (a_1, ..., a_n, a_n, ..., a_1)
reversed_part = rational_part[::-1]
period = rational_part + reversed_part
final_cf_str = f"[0; ({', '.join(map(str, period))})]"
print(f"The associated continued fraction is {final_cf_str}.")

# 4. Compute and print the quadratic equation for its value.
A, B, C = get_quadratic_equation_from_period(period)

# Format the equation string for printing
equation_str = f"{A}*x^2"
equation_str += f" + {B}*x" if B > 0 else f" - {-B}*x" if B < 0 else ""
equation_str += f" + {C}" if C > 0 else f" - {-C}" if C < 0 else ""
equation_str += " = 0"

print("\nThe value of this continued fraction, x, is the positive root of the quadratic equation:")
print(equation_str)

print("\nThe numbers in the final equation are:")
print(A)
print(B)
print(C)