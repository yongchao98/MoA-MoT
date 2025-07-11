import math
from itertools import combinations_with_replacement

def get_monotone_maps(domain_size, codomain):
    """Generates monotone maps from {0, ..., domain_size-1} to codomain."""
    maps = []
    # A map is defined by the list of its values [f(0), f(1), ..., f(domain_size-1)]
    for values in combinations_with_replacement(codomain, domain_size):
        f = list(values)
        maps.append(f)
    return maps

def solve():
    """
    Calculates the cardinality of [F*G](bullet) for n=8.
    This is equivalent to counting the number of monoid homomorphisms from A to A.
    """
    
    # The set V of exponent vectors for numbers 1..8
    # Primes are 2, 3, 5, 7. A number x is represented by (v2(x), v3(x), v5(x), v7(x)).
    V_set = {
        (0, 0, 0, 0),  # 1
        (1, 0, 0, 0),  # 2
        (0, 1, 0, 0),  # 3
        (2, 0, 0, 0),  # 4
        (0, 0, 1, 0),  # 5
        (1, 1, 0, 0),  # 6
        (0, 0, 0, 1),  # 7
        (3, 0, 0, 0)   # 8
    }
    
    inf_vec = (float('inf'),) * 4
    V_with_inf = V_set.union({inf_vec})

    # Case 1: f(1) = infinity.
    # This implies f(x) = infinity for all x in A. This is one homomorphism.
    # f(gcd(a,b))=inf. gcd(f(a),f(b))=gcd(inf,inf)=inf. f(inf)=inf.
    count = 1
    print("Counting homomorphisms f where f(1) = infinity: 1")

    # Case 2: f(1) = 1.
    # This implies h_p(0) = 0 for all primes p.
    # Let's generate the monotone functions h_p for p=2,3,5,7.
    # h_p maps v_p(x) to v_p(f(x)).
    # h_p(inf) = inf. We handle this by using float('inf').
    
    codomain_h2 = [0, 1, 2, 3, float('inf')]
    codomain_h3 = [0, 1, float('inf')]
    codomain_h5 = [0, 1, float('inf')]
    codomain_h7 = [0, 1, float('inf')]

    # Generate maps h_p where h_p(0)=0.
    # Domain for h2 is {0,1,2,3}. So we fix f(0)=0 and choose f(1),f(2),f(3).
    h2_maps = []
    for c in combinations_with_replacement(codomain_h2, 3):
        # c = (y1, y2, y3) where y1=h2(1), y2=h2(2), y3=h2(3)
        # We need to enforce monotonicity with h2(0)=0, so 0 <= y1 <= y2 <= y3
        h2_map = {0: 0, 1: c[0], 2: c[1], 3: c[2]}
        h2_maps.append(h2_map)

    # Domain for h3 is {0,1}. We fix f(0)=0 and choose f(1).
    h3_maps = []
    for c in codomain_h3:
        h3_map = {0: 0, 1: c}
        h3_maps.append(h3_map)
        
    h5_maps = h3_maps
    h7_maps = h3_maps

    num_h2 = len(h2_maps)
    num_h3 = len(h3_maps)
    num_h5 = len(h5_maps)
    num_h7 = len(h7_maps)

    print(f"Generated {num_h2} monotone maps for h_2 with h_2(0)=0.")
    print(f"Generated {num_h3} monotone maps for h_3 with h_3(0)=0.")
    print(f"Generated {num_h5} monotone maps for h_5 with h_5(0)=0.")
    print(f"Generated {num_h7} monotone maps for h_7 with h_7(0)=0.")

    valid_f_count = 0
    for h2 in h2_maps:
        for h3 in h3_maps:
            for h5 in h5_maps:
                for h7 in h7_maps:
                    is_valid_homomorphism = True
                    for v in V_set:
                        # Apply the homomorphism in the exponent space
                        v_image = (h2[v[0]], h3[v[1]], h5[v[2]], h7[v[3]])
                        
                        # Check if the image vector is in V or is the infinity vector
                        if v_image not in V_with_inf:
                            is_valid_homomorphism = False
                            break
                    
                    if is_valid_homomorphism:
                        valid_f_count += 1
            
    print(f"\nCounting homomorphisms f where f(1) = 1: {valid_f_count}")

    total_count = count + valid_f_count
    
    print(f"\nThe number of homomorphisms f with f(1)=infinity is 1.")
    print(f"The number of homomorphisms f with f(1)=1 is {valid_f_count}.")
    print(f"The total number of homomorphisms is {1} + {valid_f_count} = {total_count}.")
    print("\nThe cardinality of [F * G](bullet) is the total number of homomorphisms.")
    print(f"Final Answer: {total_count}")

solve()