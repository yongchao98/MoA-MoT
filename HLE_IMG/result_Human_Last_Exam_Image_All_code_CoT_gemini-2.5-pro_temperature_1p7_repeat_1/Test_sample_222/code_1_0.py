import math
from fractions import Fraction

def solve_contour_integral():
    """
    Calculates the imaginary part of the sum of the integrals of f(z) over contours C1 and C2.
    f(z) = z / (z - 3/2) * Gamma(z)
    """

    # Step 1 & 2: Identify singularities and calculate their residues.
    # The singularities of f(z) are simple poles at z=1.5 and at the non-positive integers z=-1, -2, -3, ...
    # The pole of Gamma(z) at z=0 is cancelled by the 'z' in the numerator, creating a removable singularity.
    
    # Residues at z = -k, k=1,2,3... : Res(f,-k) = (z/(z-1.5)) * Res(Gamma, -k)
    # Res(Gamma, -k) = (-1)^k / k!
    # So, Res(f,-k) = (-k/(-k-1.5)) * (-1)^k / k! = (k/(k+1.5)) * (-1)^k / k! = (2k/(2k+3)) * (-1)^k / k!
    
    res_minus_1 = -2/5.0      # For k=1: (2/5) * (-1)/1
    res_minus_2 = 2/7.0       # For k=2: (4/7) * 1/2
    res_minus_3 = -1/9.0      # For k=3: (6/9) * (-1)/6
    
    # Residue at z=1.5: Res(f,1.5) = z*Gamma(z)|_{z=1.5} = 1.5 * Gamma(1.5)
    # Gamma(1.5) = Gamma(1/2 + 1) = (1/2)*Gamma(1/2) = sqrt(pi)/2
    res_1_5 = 1.5 * math.sqrt(math.pi) / 2.0  # (3/2) * (sqrt(pi)/2) = 3*sqrt(pi)/4
    
    residues_str = {
        -1: "-2/5", -2: "2/7", -3: "-1/9", 1.5: "3*sqrt(pi)/4"
    }
    
    # Step 3: Determine winding numbers for each contour around each pole.
    # From visual inspection of the image:
    winding_C1 = { -1: -1, -2: 1, -3: -1, 1.5: 1 } # CW, CCW, CW, CCW loops
    winding_C2 = { -1: 1,  -2: 1, -3: 0, 1.5: 0 }  # CCW simple loop
    
    poles = [-3, -2, -1, 1.5]
    
    # Step 4: Calculate total winding numbers by summing C1 and C2.
    total_winding_numbers = {p: winding_C1[p] + winding_C2[p] for p in poles}

    # Step 5 & 6: Apply the residue theorem and find the imaginary part.
    # The sum of integrals is I = 2*pi*i * Sum[ n_total(zk) * Res(f, zk) ]
    # The term inside the sum, S = Sum[ n_total(zk) * Res(f, zk) ], is real.
    # Im(I) = 2*pi*S
    
    # Calculate S from its constituent parts.
    term_minus_3 = total_winding_numbers[-3] * res_minus_3
    term_minus_2 = total_winding_numbers[-2] * res_minus_2
    term_minus_1 = total_winding_numbers[-1] * res_minus_1
    term_1_5 = total_winding_numbers[1.5] * res_1_5
    
    # S is the sum of these terms
    S_val = term_minus_3 + term_minus_2 + term_minus_1 + term_1_5
    
    # The final imaginary part is 2 * pi * S
    imaginary_part_val = 2 * math.pi * S_val

    # Print out the derivation
    print("The sum of the integrals is I = I(C1) + I(C2). By the Residue Theorem, this is:")
    print("I = 2*pi*i * Sum[ (n(C1,zk) + n(C2,zk)) * Res(f,zk) ]")
    print("\nThe total winding numbers n_total(zk) = n(C1,zk) + n(C2,zk) for each pole zk are:")
    print(f"n_total(-3) = {winding_C1[-3]} + {winding_C2[-3]} = {total_winding_numbers[-3]}")
    print(f"n_total(-2) = {winding_C1[-2]} + {winding_C2[-2]} = {total_winding_numbers[-2]}")
    print(f"n_total(-1) = {winding_C1[-1]} + {winding_C2[-1]} = {total_winding_numbers[-1]}")
    print(f"n_total(1.5) = {winding_C1[1.5]} + {winding_C2[1.5]} = {total_winding_numbers[1.5]}")

    print("\nThe sum term S = Sum[ n_total(zk) * Res(f,zk) ] is:")
    print(f"S = ({total_winding_numbers[-3]})*({residues_str[-3]}) + ({total_winding_numbers[-2]})*({residues_str[-2]}) + ({total_winding_numbers[-1]})*({residues_str[-1]}) + ({total_winding_numbers[1.5]})*({residues_str[1.5]})")
    print("S = (1/9) + (4/7) + 0 + 3*sqrt(pi)/4")
    
    rat_sum = Fraction(1,9) + Fraction(4,7)
    print(f"S = {rat_sum.numerator}/{rat_sum.denominator} + 3*sqrt(pi)/4")

    print("\nThe imaginary part of the sum of integrals is Im(I) = 2*pi*S.")
    print("So, the final equation for the imaginary part is:")
    print(f"Im(I) = 2 * pi * ( (1/9) + (4/7) + (3*sqrt(pi)/4) )")
    
    print("\nSimplifying the expression:")
    print(f"Im(I) = 2 * pi * ( {rat_sum.numerator}/{rat_sum.denominator} + 3*sqrt(pi)/4 )")
    print(f"Im(I) = {2*rat_sum.numerator}*pi/{rat_sum.denominator} + {2*3}*pi*sqrt(pi)/4")
    print(f"Im(I) = {86}*pi/{63} + {3}*pi*sqrt(pi)/{2}")

    print(f"\nThe numerical value is: {imaginary_part_val}")
    
    return imaginary_part_val

# Run the solver and print the final answer in the required format
final_answer = solve_contour_integral()
# The instruction asks for the format <<<answer content>>> at the end
# e.g., <<<C>>>, <<<9.8>>>
# So, I'll output the numerical result rounded to one decimal place.
print(f"\n<<<{final_answer:.1f}>>>")
