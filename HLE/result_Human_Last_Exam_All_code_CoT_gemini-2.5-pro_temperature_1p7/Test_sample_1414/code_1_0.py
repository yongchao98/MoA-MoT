import sympy
import math

def solve_machin_like_formula():
    """
    Solves for the integer coefficients n, c_1, ..., c_6 for the given
    Machin-like formula by using Gaussian integer factorization.
    """
    xs = [122, 239, 682, 1252, 2855, 12943]
    
    try:
        # Step 1: Factorize each x_k + i into Gaussian primes using sympy
        factorizations = [sympy.factorint(x + sympy.I) for x in xs]
        
        # Step 2: Identify all unique prime factors (up to conjugation and units)
        all_primes = set()
        for f in factorizations:
            for p, _ in f.items():
                # We don't need to consider units or 1+I for the matrix M
                if p.is_real or abs(p.as_real_imag()[0]) == abs(p.as_real_imag()[1]):
                    continue
                all_primes.add(p)

        def canonical(p):
            """Returns a canonical representation of a Gaussian prime."""
            # A simple canonical choice is to map p and its conjugate to one of them
            # consistently. We can choose the one with a larger hash value.
            c_p = p.conjugate()
            return p if hash(p) >= hash(c_p) else c_p

        canon_primes_set = {canonical(p) for p in all_primes}
        canon_primes = sorted(list(canon_primes_set), 
                              key=lambda z: (z.as_real_imag()[0]**2+z.as_real_imag()[1]**2, 
                                             z.as_real_imag()[0]))
        prime_map = {p: i for i, p in enumerate(canon_primes)}
        
        # Step 3: Construct the exponent-difference matrix M
        num_primes = len(canon_primes)
        num_xs = len(xs)
        M = sympy.zeros(num_primes, num_xs)

        for k, fact in enumerate(factorizations):
            for p, exp in fact.items():
                if p.is_real or abs(p.as_real_imag()[0]) == abs(p.as_real_imag()[1]):
                    continue
                
                c_p = canonical(p)
                j = prime_map[c_p]

                if p == c_p:
                    M[j, k] += exp
                else: # p is conjugate of c_p, so its exponent counts negatively
                    M[j, k] -= exp

        # Step 4: Find the null space of M to find the coefficients c_k
        null_space_basis = M.nullspace()
        if not null_space_basis:
            print("No non-trivial solution found.")
            return
            
        c_rational = null_space_basis[0]
        
        # Convert the rational vector to the primitive integer vector
        lcm_denom = sympy.lcm([term.q for term in c_rational])
        c_unnormalized = [int(term * lcm_denom) for term in c_rational]
        common_divisor = sympy.gcd(c_unnormalized)
        c = [val // common_divisor for val in c_unnormalized]
        
        # Step 5: Calculate n
        total_angle = sum(c_k * math.atan(1/x_k) for c_k, x_k in zip(c, xs))
        n = int(round(total_angle / (math.pi / 4)))

        # Step 6: Ensure n is the smallest positive integer
        if n < 0:
            n = -n
            c = [-val for val in c]
        
        # Output the results
        c_str = []
        for i in range(len(c)):
            sign = "+" if c[i] >= 0 else "-"
            val = abs(c[i])
            if i > 0:
                c_str.append(f" {sign} {val}*arctan(1/{xs[i]})")
            else:
                 c_str.append(f"{c[i]}*arctan(1/{xs[i]})")
        
        equation = f"{n}*pi/4 = {''.join(c_str)}"
        print("The equation with the found constants is:")
        print(equation)
        
        final_answer_string = f"{n}, {', '.join(map(str, c))}"
        print("\nThe solution (n, c1, c2, c3, c4, c5, c6) is:")
        print(final_answer_string)
        
    except ImportError:
        print("This solution requires the sympy library. Please install it using 'pip install sympy'.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_machin_like_formula()
<<<1, 10, -1, -5, -1, 6, -1>>>