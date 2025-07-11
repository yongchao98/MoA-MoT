import math

def calculate_area():
    """
    Calculates the area of the region R based on the derived formulas.
    """
    pi = math.pi
    sqrt = math.sqrt
    asin = math.asin

    # The total area is the sum of areas from different integer pairs (n, m).
    # Area = A_00 + (A_01+A_10) + (A_02+A_20) + ... + (A_34+A_43)

    # Term for (n,m) = (0,0) -> k=0
    # Area is a quarter circle of radius 1. Formula: pi/4
    term_00 = pi / 4
    
    # Term for (n,m) = (0,1), (1,0) -> k=1
    # Formula: 2 * (Integral from 0 to 1 of sqrt(4-x^2) dx - 1)
    # Integral part evaluates to: (sqrt(3)/2 + 2*asin(1/2)) = sqrt(3)/2 + pi/3
    term_01 = 2 * (sqrt(3)/2 + pi/3 - 1)

    # Term for (n,m) = (0,2), (2,0) -> k=2
    # Formula: 2 * (Integral from 0 to 1 of sqrt(9-x^2) dx - 2)
    # Integral part evaluates to: (sqrt(8)/2 + 9/2*asin(1/3)) = sqrt(2) + 9/2*asin(1/3)
    term_02 = 2 * (sqrt(2) + (9/2)*asin(1/3) - 2)
    
    # Term for (n,m) = (0,3), (3,0) -> k=3
    # Formula: 2 * (Integral from 0 to 1 of sqrt(16-x^2) dx - 3)
    # Integral part evaluates to: (sqrt(15)/2 + 16/2*asin(1/4))
    term_03 = 2 * (sqrt(15)/2 + 8*asin(1/4) - 3)

    # Term for (n,m) = (0,4), (4,0) -> k=4
    # Formula: 2 * (Integral from 0 to 1 of sqrt(25-x^2) dx - 4)
    # Integral part evaluates to: (sqrt(24)/2 + 25/2*asin(1/5))
    term_04 = 2 * (sqrt(24)/2 + (25/2)*asin(1/5) - 4)

    # For k=5, we have two sets of pairs: ((0,5), (5,0)) and ((3,4), (4,3))
    
    # Term for (n,m) = (0,5), (5,0) -> k=5
    # Formula: 2 * (Integral from 0 to 1 of sqrt(36-x^2) dx - 5)
    # Integral part evaluates to: (sqrt(35)/2 + 36/2*asin(1/6))
    term_05 = 2 * (sqrt(35)/2 + 18*asin(1/6) - 5)
    
    # Term for (n,m) = (3,4), (4,3) -> k=5
    # Formula: 2 * (Integral from 3 to 4 of sqrt(36-x^2) dx - 4)
    # Integral part evaluates to: (2*sqrt(20)+18*asin(2/3)) - (9*sqrt(3)/2 + 3*pi)
    term_34 = 2 * (4*sqrt(5) + 18*asin(2/3) - (9*sqrt(3)/2) - 3*pi - 4)

    # Print the breakdown of the calculation
    print("The total area is the sum of the following components:")
    print(f"Contribution from (0,0): A_00 = pi/4 = {term_00:.4f}")
    print(f"Contribution from (0,1) & (1,0): A_01+A_10 = sqrt(3) + 2*pi/3 - 2 = {term_01:.4f}")
    print(f"Contribution from (0,2) & (2,0): A_02+A_20 = 2*sqrt(2) + 9*asin(1/3) - 4 = {term_02:.4f}")
    print(f"Contribution from (0,3) & (3,0): A_03+A_30 = sqrt(15) + 16*asin(1/4) - 6 = {term_03:.4f}")
    print(f"Contribution from (0,4) & (4,0): A_04+A_40 = 2*sqrt(6) + 25*asin(1/5) - 8 = {term_04:.4f}")
    print(f"Contribution from (0,5) & (5,0): A_05+A_50 = sqrt(35) + 36*asin(1/6) - 10 = {term_05:.4f}")
    print(f"Contribution from (3,4) & (4,3): A_34+A_43 = 8*sqrt(5) + 36*asin(2/3) - 9*sqrt(3) - 6*pi - 8 = {term_34:.4f}")
    
    # Total Area
    total_area = term_00 + term_01 + term_02 + term_03 + term_04 + term_05 + term_34
    
    print("\nTotal Area Calculation:")
    print(f"{term_00:.4f} + {term_01:.4f} + {term_02:.4f} + {term_03:.4f} + {term_04:.4f} + {term_05:.4f} + {term_34:.4f} = {total_area:.4f}")
    
    print(f"\nThe total area of R is: {total_area:.2f}")

calculate_area()