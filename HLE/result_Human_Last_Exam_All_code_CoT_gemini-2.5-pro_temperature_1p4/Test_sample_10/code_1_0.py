import sympy as sp

def solve():
    """
    This function calculates the derivative of T2(x) and prints it,
    highlighting the coefficients in the final equation.
    """
    x = sp.symbols('x')
    
    # Define T2(x)
    # T2(x) = ((-1 + (1 + sp.exp(x))**2) * x) / (1 + (1 + sp.exp(x))**2)
    # Simplified T2(x) = x * g(x)
    
    u = sp.exp(x)
    
    # g(x) is the fraction part
    g_x = (2*u + u**2) / (2 + 2*u + u**2)
    
    # g'(x) is the derivative of g(x)
    g_prime_x = (4*u + 4*u**2) / (2 + 2*u + u**2)**2
    
    # T2'(x) = g(x) + x * g'(x)
    
    print("The function whose derivative is not written in terms of sigmoid is T2(x).")
    print("Its derivative is T2'(x) = g(x) + x*g'(x), where:")
    
    # Printing the final equation with numbers explicitly shown
    print("g(x) = ( (2)*e^x + (1)*e^(2x) ) / ( (2) + (2)*e^x + (1)*e^(2x) )")
    print("g'(x) = ( (4)*e^x + (4)*e^(2x) ) / ( (2) + (2)*e^x + (1)*e^(2x) )^2")

solve()
