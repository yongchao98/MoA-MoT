import math

def solve_overhang():
    """
    Calculates the coefficients a, b, c for the maximal overhang of three cubes.

    The maximal overhang for three identical cubes, allowing for rotation, is
    H = 11 / (6 * sqrt(2)) block sidelengths.

    We need to express this in the form (a + sqrt(b)) / (1 + c).

    H = 11 * sqrt(2) / 12
      = sqrt(11^2 * 2) / 12
      = sqrt(242) / 12
      = (0 + sqrt(242)) / 12
      = (0 + sqrt(242)) / (1 + 11)

    This gives:
    a = 0
    b = 242
    c = 11

    The conditions are met:
    - a, b, c are non-negative integers.
    - sqrt(b) = sqrt(242) is non-integer.
    - c = 11 is minimal for the denominator 12.
    """

    # Numerator of the fraction part of the sum of harmonic series-like terms (1/2 + 1/3 + 1/1)
    # This comes from the centers of mass calculation
    # H = (1 + 1/2 + 1/3)*R_max = (11/6)*R_max
    # R_max = 1/sqrt(2)
    # H = (11/6)/sqrt(2) = 11/(6*sqrt(2))
    
    # We are asked to provide a, b, c for the expression (a + sqrt(b))/(1+c)
    # Overhang = 11 / (6 * sqrt(2)) = 11*sqrt(2)/12 = sqrt(121*2)/12 = sqrt(242)/12
    # So a=0, b=242, c=11

    a = 0
    b = 11**2 * 2
    c = 12 - 1
    
    # The maximal overhang is (0 + sqrt(242)) / (1 + 11)
    # The code prints the values for a, b, and c as requested.
    # The problem requires to "output each number in the final equation!".
    # Given the examples like "11 0 11", we print the numbers separated by spaces.
    print(f"{a} {b} {c}")

solve_overhang()