import math

def calculate_area():
    """
    Calculates the area of the region R based on the analytical solution.
    """
    total_area = 0.0
    
    print("Calculating the area of region R. The total area is a sum of contributions from several unit squares.")
    print("-----------------------------------------------------------------")

    # Case (0,0): ka=0, kb=0, K=0. Area is for a^2+b^2 < 1 in [0,1)x[0,1).
    # This is a quarter circle of radius 1.
    # Area = (1/4) * pi * 1^2
    a_0_0 = math.pi / 4.0
    total_area += a_0_0
    print(f"Contribution from square (0,0): A(0,0) = pi/4 = {a_0_0:.4f}")

    # Cases (0,k) and (k,0) for k=1..5
    # For (0,1): ka=0, kb=1, K=1. Area is for a^2+b^2 < 4 in [0,1)x[1,2).
    # Integral from 0 to 1 of (sqrt(4-a^2) - 1) da
    a_0_1 = math.sqrt(3)/2.0 + math.pi/3.0 - 1.0
    total_area += 2 * a_0_1
    print(f"Contribution from squares (0,1) and (1,0): 2 * A(0,1) = 2 * {a_0_1:.4f} = {2*a_0_1:.4f}")

    # For (0,2): ka=0, kb=2, K=2. Area is for a^2+b^2 < 9 in [0,1)x[2,3).
    # Integral from 0 to 1 of (sqrt(9-a^2) - 2) da
    a_0_2 = math.sqrt(2) + 4.5 * math.asin(1.0/3.0) - 2.0
    total_area += 2 * a_0_2
    print(f"Contribution from squares (0,2) and (2,0): 2 * A(0,2) = 2 * {a_0_2:.4f} = {2*a_0_2:.4f}")

    # For (0,3): ka=0, kb=3, K=3. Area is for a^2+b^2 < 16 in [0,1)x[3,4).
    # Integral from 0 to 1 of (sqrt(16-a^2) - 3) da
    a_0_3 = math.sqrt(15)/2.0 + 8.0 * math.asin(1.0/4.0) - 3.0
    total_area += 2 * a_0_3
    print(f"Contribution from squares (0,3) and (3,0): 2 * A(0,3) = 2 * {a_0_3:.4f} = {2*a_0_3:.4f}")

    # For (0,4): ka=0, kb=4, K=4. Area is for a^2+b^2 < 25 in [0,1)x[4,5).
    # Integral from 0 to 1 of (sqrt(25-a^2) - 4) da
    a_0_4 = math.sqrt(6) + 12.5 * math.asin(1.0/5.0) - 4.0
    total_area += 2 * a_0_4
    print(f"Contribution from squares (0,4) and (4,0): 2 * A(0,4) = 2 * {a_0_4:.4f} = {2*a_0_4:.4f}")

    # For (0,5): ka=0, kb=5, K=5. Area is for a^2+b^2 < 36 in [0,1)x[5,6).
    # Integral from 0 to 1 of (sqrt(36-a^2) - 5) da
    a_0_5 = math.sqrt(35)/2.0 + 18.0 * math.asin(1.0/6.0) - 5.0
    total_area += 2 * a_0_5
    print(f"Contribution from squares (0,5) and (5,0): 2 * A(0,5) = 2 * {a_0_5:.4f} = {2*a_0_5:.4f}")

    # Case (3,4) and (4,3): ka=3, kb=4, K=5. Area is for a^2+b^2 < 36 in [3,4)x[4,5).
    # Integral from 3 to 4 of (min(5, sqrt(36-a^2)) - 4) da
    a_3_4 = 16.0 - 2.5*math.sqrt(11) - 4.5*math.sqrt(3) - 3.0*math.pi + 18.0*math.asin(math.sqrt(11)/6.0)
    total_area += 2 * a_3_4
    print(f"Contribution from squares (3,4) and (4,3): 2 * A(3,4) = 2 * {a_3_4:.4f} = {2*a_3_4:.4f}")
    
    print("-----------------------------------------------------------------")
    print(f"The total area is the sum of all these contributions.")
    print(f"Total Area = {a_0_0:.4f} + {2*a_0_1:.4f} + {2*a_0_2:.4f} + {2*a_0_3:.4f} + {2*a_0_4:.4f} + {2*a_0_5:.4f} + {2*a_3_4:.4f}")
    print(f"Total Area = {total_area:.4f}")
    print(f"The area of R expressed to two decimals is: {total_area:.2f}")

calculate_area()