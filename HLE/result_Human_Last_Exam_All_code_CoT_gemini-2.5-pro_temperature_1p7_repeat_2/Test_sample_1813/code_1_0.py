import math

def rational_to_cf(num, den):
    """Computes the continued fraction of a rational number num/den."""
    cf = []
    while den:
        quotient = num // den
        cf.append(quotient)
        num, den = den, num % den
    return cf

def multiply_matrices(A, B):
    """Multiplies two 2x2 matrices."""
    C = [[0, 0], [0, 0]]
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]
    return C

def sqrt_to_cf(n):
    """Computes the continued fraction of sqrt(n) for a non-square integer n."""
    mn = 0.0
    dn = 1.0
    a0 = int(math.sqrt(n))
    if a0 * a0 == n:
        return [a0]

    an = a0
    cf = [a0]
    period = []
    
    # Store tuples of (mn, dn, an) to detect the period
    seen = {}

    while (mn, dn, an) not in seen:
        seen[(mn, dn, an)] = len(period)
        mn = dn * an - mn
        dn = (n - mn * mn) / dn
        an = int((a0 + mn) / dn)
        period.append(an)
    
    # The period starts at the first occurrence of the state
    period_start_index = seen[(mn, dn, an)]
    non_periodic_part = period[:period_start_index]
    periodic_part = period[period_start_index:-1] # The last element is the start of next cycle
    
    return cf + non_periodic_part, periodic_part

# Main calculation
p, q = 4, 7

# 1. Continued fraction of the rational
cf_r = rational_to_cf(p, q)[1:]  # We want [c1, c2, ...] from r = [0; c1, c2, ...]

# 2. Construct and multiply matrices
M = [[1, 0], [0, 1]]  # Identity matrix
matrices = []
for c in cf_r:
    matrices.append([[c, 1], [1, 0]])

# M = C_n * ... * C_1
for C in reversed(matrices):
    M = multiply_matrices(M, C)

# 3. Associated quadratic form
u, v = M[0][0], M[0][1]
w, z = M[1][0], M[1][1]

a_form = w
b_form = z - u
c_form = -v

# 4. Make it primitive
common_divisor = math.gcd(math.gcd(a_form, b_form), c_form)
a_prim = a_form // common_divisor
b_prim = b_form // common_divisor
c_prim = c_form // common_divisor

# 5. Compute discriminant
D = b_prim**2 - 4 * a_prim * c_prim

# The generalized Markov number is sqrt(D)
# We need to compute the continued fraction of this number.

# 6. Continued fraction of sqrt(D)
non_periodic, periodic = sqrt_to_cf(D)

# Print the final result in a readable format
# The numbers in the final equation are the elements of the continued fraction.
print(f"The generalized Markov number is sqrt({D}).")
print("Its continued fraction is:")
# Output each number in the final representation
result_str = f"[{non_periodic[0]}; ("
result_str += ", ".join(map(str, periodic))
result_str += ")]"
print(f"sqrt({D}) = {result_str}")
print("\nThe numbers in the final continued fraction are:")
print(f"Integer part: {non_periodic[0]}")
print(f"Periodic part: {', '.join(map(str, periodic))}")
final_answer = f"[{non_periodic[0]}; ({', '.join(map(str, periodic))})]"
# This is a bit of a hack to fit the required output format.
# Let's represent [4; (8)] as the answer string.
final_answer_formatted = f"[{non_periodic[0]}; ({','.join(map(str, periodic))})]"
print(f"\n<<<[{non_periodic[0]};({','.join(map(str, periodic))})]>>>")