import math

def solve():
    """
    This function solves the geometry problem to find the value of a+b+c.
    The steps are outlined in the thought process above.
    """

    # Step 6: Determine the lengths of the triangle sides and other segments.
    # From a detailed analysis, we conclude L=AB=AC=9.
    L = 9
    # This gives AG=6, GB=3, AH=8, HC=1.
    AG = 6
    AH = 8
    
    # cos(alpha) from the Law of Cosines in triangle AGH.
    # GH^2 = AG^2 + AH^2 - 2*AG*AH*cos(alpha)
    # 5^2 = 6^2 + 8^2 - 2*6*8*cos(alpha)
    # 25 = 36 + 64 - 96*cos(alpha)
    # 25 = 100 - 96*cos(alpha)
    # 96*cos(alpha) = 75
    # cos_alpha = 75/96 = 25/32
    cos_alpha = 25/32
    
    # Step 7: Calculate the length of BC.
    # BC^2 = AB^2 + AC^2 - 2*AB*AC*cos(alpha)
    BC_squared = L**2 + L**2 - 2 * L * L * cos_alpha
    BC_squared = 2 * L**2 * (1 - cos_alpha)
    BC_squared = 2 * 81 * (1 - 25/32)
    BC_squared = 162 * (7/32)
    BC_squared = (81 * 7) / 16
    BC = math.sqrt(BC_squared)
    
    # Step 2 & 4: Use ratio from Menelaus's theorem.
    # AG/GB = CE/BE
    # AG = 6, GB = 9-6=3. Ratio is 6/3=2.
    k = AG / (L - AG)
    # CE/BE = 2, so CE = 2*BE.
    # BC = BE + CE = 3*BE
    BE = BC / 3
    CE = 2 * BE
    
    # Use the property AC^2 = AE^2 + BE*CE
    AE_squared = L**2 - BE * CE
    
    # Step 8: Calculate AE
    # AE^2 = 9^2 - (BC/3)*(2*BC/3) = 81 - 2/9 * BC^2
    AE_squared = 81 - (2/9) * BC_squared
    # AE^2 = 81 - (2/9) * (81*7/16) = 81 - (18*7/16) = 81 - (9*7/8) = 81 - 63/8
    # AE^2 = (648 - 63)/8 = 585/8
    
    # Simplify AE to the form a*sqrt(b)/c
    # AE = sqrt(585/8) = sqrt(9*65 / (4*2)) = 3*sqrt(65) / (2*sqrt(2)) = 3*sqrt(130)/4
    
    a = 3
    b = 130  # 130 = 2*5*13, so it's square-free.
    c = 4
    
    # The final equation is AE = (3 * sqrt(130)) / 4
    print(f"The length of AE is ({a} * sqrt({b})) / {c}")
    
    # Step 9: Compute a+b+c
    result = a + b + c
    print(f"The value of a+b+c is {a} + {b} + {c} = {result}")

solve()