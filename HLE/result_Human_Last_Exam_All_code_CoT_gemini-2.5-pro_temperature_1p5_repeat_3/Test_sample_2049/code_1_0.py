import math

def get_divisors(n):
    """Helper function to get all divisors of a number n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return list(divs)

def solve():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.
    This is equivalent to counting the number of etale algebras of degree 5
    over Q unramified outside {2, infinity}.
    """
    # N[k] is the number of number fields of degree k unramified outside {2, infinity}.
    # These values are from the LMFDB database.
    N = {
        1: 1,  # Q
        2: 3,  # Q(i), Q(sqrt(2)), Q(sqrt(-2))
        3: 0,
        4: 10,
        5: 6,
    }
    
    # A[n] will store the number of etale algebras of degree n.
    # We initialize with A[0] = 1 (the trivial algebra of degree 0).
    A = {0: 1}
    
    # G[k] = sum_{d|k} d * N(d)
    G = {}
    
    # We want to compute A(5). We compute A(n) for n from 1 to 5.
    target_degree = 5
    for n in range(1, target_degree + 1):
        sum_val = 0
        for k in range(1, n + 1):
            # Calculate G[k] if not already computed
            if k not in G:
                G_k = 0
                divisors = get_divisors(k)
                for d in divisors:
                    G_k += d * N.get(d, 0)
                G[k] = G_k
            
            sum_val += G[k] * A[n-k]
        
        A[n] = int(sum_val / n)

    # Print the final calculation for A(5)
    # 5 * A(5) = G(1)A(4) + G(2)A(3) + G(3)A(2) + G(4)A(1) + G(5)A(0)
    
    terms = []
    total_sum = 0
    for k in range(1, target_degree + 1):
        term_val = G[k] * A[target_degree-k]
        terms.append(f"{G[k]} * {A[target_degree-k]}")
        total_sum += term_val
    
    equation = f"{target_degree} * A({target_degree}) = {' + '.join(terms)}"
    print(f"The number of etale algebras A(n) is computed recursively.")
    for i in range(1, target_degree):
        print(f"A({i}) = {A[i]}")

    print("\nThe final calculation for A(5) is based on the formula:")
    print(f"5 * A(5) = G(1)*A(4) + G(2)*A(3) + G(3)*A(2) + G(4)*A(1) + G(5)*A(0)")
    
    print("\nSubstituting the computed values:")
    g_values_str = ', '.join([f'G({i})={G[i]}' for i in range(1,target_degree+1)])
    print(f"where {g_values_str}")
    
    final_sum_str = equation.replace(f"{target_degree} * A({target_degree}) = ", "")
    print(f"{target_degree} * A({target_degree}) = {final_sum_str} = {total_sum}")

    final_result = A[target_degree]
    print(f"A({target_degree}) = {total_sum} / {target_degree} = {final_result}")
    
    print("\nThus, the number of isomorphism classes of del Pezzo surfaces of degree 5 over Q with good reduction outside 2 is:")
    print(final_result)

solve()