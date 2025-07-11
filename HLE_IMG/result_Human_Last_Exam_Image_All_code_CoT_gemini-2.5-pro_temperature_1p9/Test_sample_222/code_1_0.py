import math
from fractions import Fraction

def calculate_imaginary_part_of_integrals():
    """
    This function calculates the imaginary part of the sum of the integrals of
    f(z) = (z / (z - 1.5)) * Gamma(z) over the contours C1 and C2.
    """

    # Step 1: Calculate residues at the relevant poles.
    # The singularity at z=0 is removable. Poles are at z=1.5, -1, -2, -3,...

    # Residue at z = 1.5
    # Res(f, 1.5) = 1.5 * Gamma(1.5) = 1.5 * (sqrt(pi)/2) = (3/4)*sqrt(pi)
    res_1_5 = 1.5 * math.sqrt(math.pi) / 2.0

    # Residue at z = -n is (2n/(2n+3)) * ((-1)^n / n!)
    # Using fractions for exact rational arithmetic.
    # For n=1, pole at z=-1
    res_neg_1 = Fraction(2 * 1, 2 * 1 + 3) * Fraction((-1)**1, math.factorial(1)) # -2/5
    # For n=2, pole at z=-2
    res_neg_2 = Fraction(2 * 2, 2 * 2 + 3) * Fraction((-1)**2, math.factorial(2)) # 2/7
    # For n=3, pole at z=-3
    res_neg_3 = Fraction(2 * 3, 2 * 3 + 3) * Fraction((-1)**3, math.factorial(3)) # -1/9

    residues = {
        1.5: res_1_5,
        -1: res_neg_1,
        -2: res_neg_2,
        -3: res_neg_3
    }

    # Step 2: Determine winding numbers from the image.
    # For C1 (blue): Arrow on right loop indicates clockwise (CW) motion.
    # Assuming alternating lobe orientations:
    # Loop at 1.5: CW -> Ind = -1
    # Loop at -1: CW -> Ind = -1
    # Loop at -2: Counter-clockwise (CCW) -> Ind = +1
    # Loop at -3: CW -> Ind = -1
    winding_numbers_C1 = {1.5: -1, -1: -1, -2: 1, -3: -1}

    # For C2 (orange): Simple CCW loop enclosing z=-1 and z=-2.
    # Ind = +1 for both.
    winding_numbers_C2 = {-1: 1, -2: 1}

    # Step 3: Calculate the sum S = Sum(Res * Ind) for both contours.
    # The integral value is 2*pi*i*S, so the imaginary part is 2*pi*S.
    
    # S1 for contour C1
    S1_rational = Fraction(0)
    S1_irrational_coeff = 0.0 # for the coefficient of sqrt(pi)

    for pole, ind in winding_numbers_C1.items():
        res = residues[pole]
        if isinstance(res, Fraction):
            S1_rational += res * ind
        else: # This is the irrational part for res_1_5
            S1_irrational_coeff += (1.5 / 2.0) * ind # The residue is (1.5/2)*sqrt(pi)

    # S2 for contour C2
    S2_rational = Fraction(0)
    for pole, ind in winding_numbers_C2.items():
        res = residues[pole]
        S2_rational += res * ind

    # Total S
    S_rational = S1_rational + S2_rational
    S_irrational_coeff = S1_irrational_coeff
    
    # S = S_rational + S_irrational_coeff * sqrt(pi)
    # The final imaginary part is 2*pi*S
    # Im(I) = 2*pi * S_rational + 2*pi * S_irrational_coeff * sqrt(pi)

    # Get numbers for the final equation printout
    term1_num = 2 * S_rational.numerator
    term1_den = S_rational.denominator

    # S_irrational_coeff = (1.5/2)*(-1) = -0.75 = -3/4
    # Second term: 2*pi * (-3/4)*sqrt(pi) = - (3/2)*pi*sqrt(pi)
    term2_num = 3
    term2_den = 2

    final_expression = f"({term1_num} * pi) / {term1_den} - ({term2_num} * pi * sqrt(pi)) / {term2_den}"
    final_value = (term1_num * math.pi) / term1_den - (term2_num * math.pi * math.sqrt(math.pi)) / term2_den

    print("The imaginary part of the sum of the integrals is given by the expression:")
    print(final_expression)
    print(f"\nThis evaluates to approximately:")
    print(final_value)


calculate_imaginary_part_of_integrals()
result = (86 * math.pi) / 63 - (3 * math.pi * math.sqrt(math.pi)) / 2
<<< (86 * pi) / 63 - (3 * pi * sqrt(pi)) / 2 >>>