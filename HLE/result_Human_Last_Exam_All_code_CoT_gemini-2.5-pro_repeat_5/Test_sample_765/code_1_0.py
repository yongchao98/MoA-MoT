import math

def F(x, R):
    """
    Calculates the integral of sqrt(R^2 - t^2) from 0 to x.
    This represents the area under the circular arc y=sqrt(R^2-t^2) from 0 to x.
    """
    if R <= 0:
        return 0
    # Clamp x to the valid range [-R, R] for asin
    if x >= R:
        x = R
    elif x <= -R:
        x = -R
    
    return (x * math.sqrt(R**2 - x**2) / 2) + (R**2 * math.asin(x / R) / 2)

def disk_area_in_origin_rect(R, x_max, y_max):
    """
    Calculates the area of the disk x^2+y^2 < R^2 in the rectangle [0, x_max] x [0, y_max].
    Assumes x_max, y_max >= 0.
    """
    if R <= 0 or x_max <= 0 or y_max <= 0:
        return 0
    
    if R**2 >= x_max**2 + y_max**2:
        return x_max * y_max
    
    # Intersection of y=sqrt(R^2-x^2) and y=y_max
    if R > y_max:
        x_s = math.sqrt(R**2 - y_max**2)
        if x_s >= x_max:
            return x_max * y_max
        else:
            return x_s * y_max + (F(x_max, R) - F(x_s, R))
    else: # R <= y_max
        if R <= x_max:
            return F(R, R) # Quarter disk area
        else: # R > x_max
            return F(x_max, R)

def disk_area_in_rect(R, x1, x2, y1, y2):
    """
    Calculates the area of the disk x^2+y^2 < R^2 in the rectangle [x1, x2] x [y1, y2].
    Uses the principle of inclusion-exclusion. Assumes rectangle is in the first quadrant.
    """
    return (disk_area_in_origin_rect(R, x2, y2) - 
            disk_area_in_origin_rect(R, x1, y2) -
            disk_area_in_origin_rect(R, x2, y1) +
            disk_area_in_origin_rect(R, x1, y1))

def get_area_contribution(na, nb):
    """
    Calculates the area contribution from the unit square starting at (na, nb).
    """
    K_sq = na**2 + nb**2
    K_float = math.sqrt(K_sq)
    
    if abs(K_float - round(K_float)) > 1e-9:
        return 0
    
    K = int(round(K_float))
    
    R_inner = K
    R_outer = K + 1
    
    x1, y1 = na, nb
    x2, y2 = na + 1, nb + 1
    
    area_outer = disk_area_in_rect(R_outer, x1, x2, y1, y2)
    area_inner = disk_area_in_rect(R_inner, x1, x2, y1, y2)
    
    return area_outer - area_inner

def main():
    """
    Main function to calculate and print the total area.
    """
    # Pairs (na, nb) with na <= nb where sqrt(na^2+nb^2) is an integer
    # for na, nb in [0, 5].
    valid_pairs_unique = [
        (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (3, 4)
    ]

    total_area = 0
    area_contributions = {}

    print("Calculating area contributions for each valid unit square region:")
    
    for na, nb in valid_pairs_unique:
        area = get_area_contribution(na, nb)
        area_contributions[(na, nb)] = area
        
        if na == nb:
            total_area += area
        else:
            # Add contributions for both (na, nb) and (nb, na) by symmetry
            total_area += 2 * area

    # Print the equation
    equation_parts = []
    print("\nFinal Area Calculation:")
    for (na, nb), area in area_contributions.items():
        if na == nb:
            print(f"Area for ({na},{nb}): {area:.4f}")
            equation_parts.append(f"{area:.4f}")
        else:
            print(f"Area for ({na},{nb}) and ({nb},{na}): 2 * {area:.4f} = {2*area:.4f}")
            equation_parts.append(f"2 * {area:.4f}")

    print("\nTotal Area = " + " + ".join(equation_parts))
    print(f"Total Area = {total_area:.2f}")

if __name__ == "__main__":
    main()