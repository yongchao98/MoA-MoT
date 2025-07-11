import sympy

def get_valuation(n, p):
    """Calculates the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    count = 0
    # sympy.Integer(n) handles potential large numbers from discriminant calculation
    n_abs = abs(sympy.Integer(n)) 
    while n_abs > 0 and n_abs % p == 0:
        count += 1
        n_abs //= p
    return count

def solve():
    """
    Calculates the thickness of the double point of the stable reduction of the given curve.
    """
    x = sympy.Symbol('x')
    # The polynomial from the curve equation z^2 = f(x)
    f = 2 * x**5 + 2 * x**3 + 1
    
    # 1. Calculate the genus of the curve
    degree = sympy.degree(f, gen=x)
    genus = (degree - 1) // 2
    
    # 2. Calculate the resultant of f and its derivative f'
    f_prime = sympy.diff(f, x)
    # The resultant is often called discriminant in software packages
    res = sympy.resultant(f, f_prime)
    
    # The discriminant in arithmetic geometry is Res(f, f') / c_d
    # where c_d is the leading coefficient.
    c_d = sympy.LC(f, x)
    discriminant = res / c_d
    
    # 3. Compute the 2-adic valuation of the discriminant
    v2_discriminant = get_valuation(discriminant, 2)
    
    # 4. Apply the formula for the thickness
    # This is a formula that holds in certain cases.
    thickness = v2_discriminant // (2 * genus)

    print(f"The equation of the curve is z^2 = {f}")
    print(f"The genus of the curve is g = ({degree} - 1) / 2 = {genus}")
    print(f"The polynomial is f(x) = {f}")
    print(f"The derivative is f'(x) = {f_prime}")
    print(f"The resultant Res(f, f') is {res}")
    print(f"The leading coefficient is {c_d}")
    print(f"The discriminant is Delta_f = {res} / {c_d} = {discriminant}")
    print(f"The 2-adic valuation of the discriminant is v2(Delta_f) = {v2_discriminant}")
    print(f"The thickness is calculated as v2(Delta_f) / (2 * g)")
    print(f"{v2_discriminant} / (2 * {genus}) = {thickness}")

solve()