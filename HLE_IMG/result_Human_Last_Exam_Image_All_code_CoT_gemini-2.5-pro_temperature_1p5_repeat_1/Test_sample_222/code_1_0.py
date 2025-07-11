import math
from fractions import Fraction

def calculate_residues():
    """Calculates the residues of the function f(z) at its poles."""
    
    # Pole at z = 1.5
    # Res(f, 1.5) = lim_{z->1.5} (z-1.5) * z/(z-1.5) * Gamma(z) = 1.5 * Gamma(1.5)
    # Gamma(1.5) = Gamma(3/2) = (1/2) * Gamma(1/2) = sqrt(pi)/2
    res_1_5 = 1.5 * math.sqrt(math.pi) / 2
    
    # Poles at z = -n for n = 1, 2, 3, ...
    # Res(f, -n) = lim_{z->-n} (z+n) * z/(z-1.5) * Gamma(z)
    #            = z/(z-1.5) |_{z=-n} * Res(Gamma, -n)
    #            = -n/(-n-1.5) * (-1)^n / n! = n/(n+1.5) * (-1)^n / n!
    
    # Pole at z = -1 (n=1)
    res_neg_1 = Fraction(1, 1 + 1.5) * Fraction((-1)**1, math.factorial(1))
    res_neg_1 = Fraction(1, 2.5) * -1
    
    # Pole at z = -2 (n=2)
    res_neg_2 = Fraction(2, 2 + 1.5) * Fraction((-1)**2, math.factorial(2))
    res_neg_2 = Fraction(2, 3.5) * Fraction(1, 2)

    # Pole at z = -3 (n=3)
    res_neg_3 = Fraction(3, 3 + 1.5) * Fraction((-1)**3, math.factorial(3))
    res_neg_3 = Fraction(3, 4.5) * Fraction(-1, 6)
    
    residues = {
        1.5: res_1_5,
        -1: res_neg_1,
        -2: res_neg_2,
        -3: res_neg_3,
    }
    return residues

def main():
    """Main function to calculate the imaginary part of the sum of integrals."""
    print("Step 1: Analyzing the function and its poles.")
    print("f(z) = z/(z-3/2) * Gamma(z)")
    print("The poles are at z = 3/2 = 1.5, and at the non-positive integers z = -1, -2, -3, ... (z=0 is a removable singularity).")
    
    poles = [1.5, -1, -2, -3]
    residues = calculate_residues()
    
    print("\nStep 2: Calculating the residues at the poles.")
    print(f"Res(f, 1.5) = 1.5 * Gamma(1.5) = (3/2) * (sqrt(pi)/2) = 3*sqrt(pi)/4 ≈ {residues[1.5]:.4f}")
    print(f"Res(f, -1) = 1/(1+1.5) * (-1)^1/1! = -1/2.5 = {residues[-1]}")
    print(f"Res(f, -2) = 2/(2+1.5) * (-1)^2/2! = 1/3.5 = {residues[-2]}")
    print(f"Res(f, -3) = 3/(3+1.5) * (-1)^3/3! = -1/9 = {residues[-3]}")

    print("\nStep 3: Determining winding numbers from the contours.")
    print("Based on the arrows and standard conventions for such diagrams:")
    # From analyzing the diagram:
    # C1 is a chain of alternating CW/CCW loops. Arrows fix the directions.
    # C2 is a simple CW loop.
    ind_c1 = {1.5: -1, -1: -1, -2: 1, -3: -1}
    ind_c2 = {1.5: 0, -1: -1, -2: -1, -3: 0}
    
    print("Winding numbers for C1:")
    for p, i in ind_c1.items(): print(f"  Ind_C1({p}) = {i}")
    print("Winding numbers for C2:")
    for p, i in ind_c2.items(): print(f"  Ind_C2({p}) = {i}")
        
    print("\nStep 4: Calculating total winding numbers for the sum of integrals I = I(C1) + I(C2).")
    total_ind = {p: ind_c1.get(p, 0) + ind_c2.get(p, 0) for p in poles}
    for p in poles:
        print(f"  Ind_total({p}) = Ind_C1({p}) + Ind_C2({p}) = {ind_c1[p]} + {ind_c2[p]} = {total_ind[p]}")

    print("\nStep 5: Applying the Residue Theorem.")
    print("Sum of integrals = 2*pi*i * Sum[ Res(f, z_k) * Ind_total(z_k) ]")
    
    sum_of_res_terms = 0
    sum_of_res_terms_symbolic = []
    
    res_1_5_term = residues[1.5] * total_ind[1.5]
    res_neg_1_term = residues[-1] * total_ind[-1]
    res_neg_2_term = residues[-2] * total_ind[-2]
    res_neg_3_term = residues[-3] * total_ind[-3]
    
    sum_of_res_terms = res_1_5_term + float(res_neg_1_term) + float(res_neg_2_term) + float(res_neg_3_term)

    print("Sum of weighted residues:")
    print(f"  Res(f, 1.5) * Ind_total(1.5) = (3*sqrt(pi)/4) * ({total_ind[1.5]})")
    print(f"  Res(f, -1) * Ind_total(-1)  = ({residues[-1]}) * ({total_ind[-1]}) = {res_neg_1_term}")
    print(f"  Res(f, -2) * Ind_total(-2)  = ({residues[-2]}) * ({total_ind[-2]}) = {res_neg_2_term}")
    print(f"  Res(f, -3) * Ind_total(-3)  = ({residues[-3]}) * ({total_ind[-3]}) = {res_neg_3_term}")
    
    print("\nThe sum of integrals is purely imaginary as all residues are real.")
    print("The imaginary part is 2*pi * (sum of weighted residues).")

    rational_part = res_neg_1_term + res_neg_2_term + res_neg_3_term
    
    print("\nFinal calculation:")
    print(f"Imaginary Part = 2 * pi * [ (3*sqrt(pi)/4) * ({total_ind[1.5]}) + ({residues[-1]}) * ({total_ind[-1]}) + ({residues[-2]}) * ({total_ind[-2]}) + ({residues[-3]}) * ({total_ind[-3]}) ]")
    print(f"Imaginary Part = 2 * pi * [ -3*sqrt(pi)/4 + {res_neg_1_term} + {res_neg_2_term} + {res_neg_3_term} ]")
    print(f"Imaginary Part = 2 * pi * [ -3*sqrt(pi)/4 + {rational_part.numerator}/{rational_part.denominator} ]")

    final_value = 2 * math.pi * sum_of_res_terms
    print(f"Imaginary Part ≈ 2 * {math.pi:.4f} * ( {sum_of_res_terms:.4f} ) ≈ {final_value:.4f}")
    print(f"<<<{final_value:.4f}>>>")

if __name__ == "__main__":
    main()
