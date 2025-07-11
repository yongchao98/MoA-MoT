import math
from fractions import Fraction
from itertools import combinations

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for the Calabi-Yau link.
    """
    # Step 1: Define weights and calculate the degree d
    w = [22, 29, 49, 50, 75]
    d = sum(w)

    print("Step 1: The weights are ({}, {}, {}, {}, {}) and the degree d = sum(w) = {}.".format(*w, d))

    # Step 2: Calculate components for the main term of the formula
    prod_w = math.prod(w)
    
    # Sum of pairwise products of weights: sum_{i<j} w_i*w_j
    sum_ww = 0
    for i, j in combinations(w, 2):
        sum_ww += i * j

    main_term_num = d * sum_ww
    main_term = Fraction(main_term_num, prod_w)
    
    print(f"Step 2: Calculated main term component.")
    print(f"  Sum of pairwise products of weights = {sum_ww}")
    print(f"  Product of all weights = {prod_w}")
    print(f"  Main term = ({d} * {sum_ww}) / {prod_w} = {main_term}")
    
    # Step 3: Calculate the orbifold correction terms
    total_correction = Fraction(0)
    print("\nStep 3: Calculating orbifold correction terms.")
    
    for i, j in combinations(range(len(w)), 2):
        wi, wj = w[i], w[j]
        g = math.gcd(wi, wj)
        
        if g > 1:
            print(f"  Found non-coprime pair (w{i+1}, w{j+1}) = ({wi}, {wj}) with gcd = {g}.")
            
            # Product of other weights
            prod_w_other = 1
            for k in range(len(w)):
                if k != i and k != j:
                    prod_w_other *= w[k]
            
            # Correction for this pair
            d_ij = Fraction(d, prod_w_other)
            correction = Fraction(g - 1, g) * d_ij
            
            print(f"    Product of other weights = {prod_w_other}")
            print(f"    Correction = ({g}-1)/{g} * ({d}/{prod_w_other}) = {correction}")
            
            total_correction += correction

    print(f"\nTotal correction term = {total_correction}")

    # Step 4: Calculate c2(X).H and the final invariant
    c2H = main_term - total_correction
    print(f"\nStep 4: Combining terms.")
    print(f"  c2(X).H = Main Term - Total Correction = {main_term} - {total_correction} = {c2H}")
    
    invariant = c2H / 12
    print(f"  Invariant = c2(X).H / 12 = {c2H} / 12 = {invariant}")

    # Step 5: Final Equation Output
    print("\nFinal Equation:")
    print(f"Invariant = ( (({d} * {sum_ww}) / {prod_w}) - {total_correction} ) / 12")
    print(f"Invariant = ( {main_term} - {total_correction} ) / 12")
    print(f"Invariant = ( {c2H} ) / 12")
    print(f"Invariant = {invariant}")
    
    return invariant

final_invariant = solve_crawley_nordstrom()
<<<{}>>>.format(final_invariant)