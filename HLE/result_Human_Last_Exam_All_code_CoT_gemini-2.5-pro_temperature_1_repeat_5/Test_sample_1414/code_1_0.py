import sympy
import math

def solve_machin_formula():
    """
    Solves for the coefficients of the given Machin-like formula.
    """
    
    # The plan is outlined in the text above.
    # Here we implement steps 4 and 5, using the pre-calculated factorizations of x_k + i.

    # Hardcoded factorizations of x_k + i for x_k in [122, 239, 682, 1252, 2855, 12943]
    # Format: { 'prime_name': exponent, ... }
    # '_conj' means the conjugate of the prime, 'p2' is 1+i, 'unit' is the complex unit.
    # Primes are identified by their norm's root, e.g., p5 for prime over 5 (2+i).
    factorizations = [
        {'p5': 1, 'p13_conj': 1, 'p229': 1, 'unit': 1},           # 122+i
        {'p2': 1, 'p13_conj': 4, 'unit': 1j},                     # 239+i
        {'p5': 3, 'p61_conj': 2, 'unit': 1},                      # 682+i
        {'p5': 1, 'p37_conj': 2, 'p229_conj': 1, 'unit': 1},      # 1252+i
        {'p2': 1, 'p13': 1, 'p37': 2, 'p229_conj': 1, 'unit': -1j},# 2855+i
        {'p2': 1, 'p5': 4, 'p13_conj': 3, 'p61_conj': 1, 'unit': -1j} # 12943+i
    ]
    primes = ['p5', 'p13', 'p37', 'p61', 'p229']
    xs = [122, 239, 682, 1252, 2855, 12943]

    # Build the matrix for the linear system sum(c_j * d_j) = 0
    # where d_j(p) = exponent of p - exponent of p_conj in factorization of x_j+i.
    matrix_rows = []
    for p in primes:
        row = []
        for f in factorizations:
            exp_p = f.get(p, 0)
            exp_p_conj = f.get(p + '_conj', 0)
            row.append(exp_p - exp_p_conj)
        matrix_rows.append(row)

    A = sympy.Matrix(matrix_rows)
    
    # Find the null space (kernel) of the matrix A to find the coefficients c_j.
    null_space = A.nullspace()

    # The null space should be one-dimensional. We get a basis vector.
    basis_vector = null_space[0]
    
    # Find the least common multiple of the denominators to get an integer vector.
    lcm = sympy.lcm([f.q for f in basis_vector])
    c_unscaled = [int(val * lcm) for val in basis_vector]

    # The smallest integer vector might have a common divisor.
    common_divisor = math.gcd(*c_unscaled)
    c_scaled = [x // common_divisor for x in c_unscaled]

    # The solution is unique, which implies a specific scaling and sign.
    # We choose the sign such that the first non-zero coefficient is positive.
    first_nonzero_idx = next((i for i, x in enumerate(c_scaled) if x), None)
    if first_nonzero_idx is not None and c_scaled[first_nonzero_idx] < 0:
        c = [-x for x in c_scaled]
    else:
        c = c_scaled
    
    # Now, determine n.
    # First, calculate n_c, the total exponent of the prime (1+i).
    v_1plus_i = [f.get('p2', 0) for f in factorizations]
    n_c = sum(c_i * v_i for c_i, v_i in zip(c, v_1plus_i))

    # Second, calculate the total complex unit of the product.
    total_unit = complex(1, 0)
    units = [f.get('unit', 1) for f in factorizations]
    for i in range(len(c)):
        total_unit *= units[i]**c[i]
        
    # The argument of the product is arg(unit) + n_c * pi/4.
    # n*pi/4 = arg(total_unit) + n_c*pi/4 + 2*k*pi for some integer k.
    # This simplifies to n = arg_unit_in_quarter_pi_turns + n_c + 8k.
    arg_unit_rad = sympy.arg(total_unit)
    arg_unit_turns = arg_unit_rad / (sympy.pi / 4)
    n_base = int(round(float(arg_unit_turns))) + n_c

    # Find the smallest positive n.
    # n = n_base (mod 8)
    n = n_base % 8
    if n == 0: # n > 0
        n = 8

    # Print the full equation as requested.
    equation_parts = []
    for i in range(len(c)):
        if c[i] == 0:
            continue
        # Use a compact form for the printout.
        if c[i] == 1:
            term = f"arctan(1/{xs[i]})"
        elif c[i] == -1:
            term = f"- arctan(1/{xs[i]})"
        elif c[i] > 0:
            term = f"+ {c[i]}*arctan(1/{xs[i]})"
        else:
            term = f"- {-c[i]}*arctan(1/{xs[i]})"
        equation_parts.append(term)
    
    # Polish the equation string
    equation = " ".join(equation_parts)
    if equation.startswith("+ "):
        equation = equation[2:]
    
    print(f"The derived identity is:")
    print(f"{n} * pi/4 = {equation}\n")

    # Print the final answer in the requested format.
    print("The constants n, c1, c2, c3, c4, c5, c6 are:")
    coeffs_str = ", ".join(map(str, c))
    print(f"{n}, {coeffs_str}")

solve_machin_formula()