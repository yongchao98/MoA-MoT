from sympy import series, Symbol, Rational

def sigma(k, n):
    """
    Computes the sum of the k-th powers of the divisors of n.
    Handle cases where n is not a positive integer, where the sum is 0.
    """
    if not isinstance(n, int) or n <= 0:
        return 0
    
    divs = []
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            divs.append(i)
            if i*i != n:
                divs.append(n//i)
    
    return sum(d**k for d in divs)

def eisenstein_q_expansion(k, n_terms):
    """
    Computes the q-expansion of the normalized Eisenstein series E_k.
    E_k(z) = 1 - (2k/B_k) * sum_{n=1 to inf} sigma_{k-1}(n) * q^n
    The factor for k=4 is 240.
    """
    q = Symbol('q')
    expansion = 1
    
    if k == 4:
        factor = 240
        sigma_k = 3
    else:
        # This implementation only requires E_4
        raise NotImplementedError("Only k=4 is implemented")
        
    for n in range(1, n_terms + 1):
        expansion += factor * sigma(sigma_k, n) * q**n
    return expansion

def solve():
    """
    Solves the problem by calculating the coefficients of the unique cusp form.
    """
    # Number of terms to compute in the q-expansions
    n_terms = 10
    q = Symbol('q')

    # 1. Get the q-expansion of E_4(z)
    E4_z = eisenstein_q_expansion(4, n_terms)

    # 2. Get the expansion of F(z) = E_4(2z)
    E4_2z = E4_z.subs(q, q**2)

    # 3. Define the form h = (E_4(z) - E_4(2z))^2, which spans the unique subspace
    #    where a_1 = 0.
    delta_2 = series(E4_z - E4_2z, q, 0, n_terms + 1)
    h = series(delta_2**2, q, 0, n_terms + 1)

    # 4. Find the first non-zero coefficient to normalize h
    first_term_power = 0
    first_term_coeff = 0
    for i in range(n_terms + 1):
        coeff = h.coeff(q, i)
        if coeff != 0:
            first_term_power = i
            first_term_coeff = coeff
            break
            
    # 5. Normalize h to get f
    f = series(h / first_term_coeff, q, 0, n_terms + 1)
    
    # 6. Get the first three non-zero coefficients of f
    coeffs = []
    current_power = first_term_power
    while len(coeffs) < 3:
        c = f.coeff(q, current_power)
        if c != 0:
            coeffs.append(c)
        current_power += 1
        
    c1, c2, c3 = coeffs
    
    # 7. Calculate the sum
    total_sum = sum(coeffs)
    
    # Print the result in the requested format
    print(f"The unique normalized cusp form is f(z) = {f.removeO()}.")
    print("Its first three non-zero coefficients are:")
    print(f"a_{first_term_power} = {c1}")
    print(f"a_{first_term_power+1} = {c2}")
    print(f"a_{first_term_power+2} = {c3}")
    print(f"The sum is: {c1} + {c2} + {c3} = {total_sum}")

solve()
<<<1031/15>>>