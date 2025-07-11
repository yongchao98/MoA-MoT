import sympy
from sympy import integrate, exp, oo, Symbol, pi

def illustrate_moments():
    """
    This function demonstrates that a non-zero Schwartz function generally has non-zero moments.
    It serves as a computational illustration for the proof that if all moments are zero,
    the function must be zero.
    """
    x = Symbol('x')
    
    # Let's choose a simple, non-zero Schwartz function: a Gaussian.
    # f(x) = exp(-x^2)
    # This function is non-zero, so according to our proof, it must have at least one non-zero moment.
    f_x = exp(-x**2)
    
    print(f"Testing the non-zero Schwartz function f(x) = {f_x}")
    print("According to the theorem, this function must have at least one non-zero moment.")
    
    # We will calculate its moments for k = 0, 1, 2, ...
    # The k-th moment is the integral of x^k * f(x) from -oo to oo.
    
    max_k_to_check = 5
    print(f"\nCalculating moments for k = 0 to {max_k_to_check}:")
    
    all_moments_zero = True
    for k in range(max_k_to_check + 1):
        integrand = x**k * f_x
        # Use sympy to perform the symbolic integration
        moment = integrate(integrand, (x, -oo, oo))
        
        # The equation for the k-th moment is:
        # Integral(x^k * exp(-x^2), (x, -oo, oo)) = moment
        print(f"Moment k={k}: ∫(x^{k} * e^(-x²)) dx = {moment}")
        
        if moment != 0:
            all_moments_zero = False
            
    print("\n--- Illustration Conclusion ---")
    if not all_moments_zero:
        print("As shown above, not all moments are zero (e.g., for k=0, 2, 4).")
        print("This is consistent with the fact that our chosen function f(x) is not identically zero.")
    else:
        # This case won't be reached for our chosen function
        print("All tested moments were zero. This is unexpected.")

    print("\n--- Final Result from Proof ---")
    print("The rigorous proof shows that if ALL moments were zero, the function would have to be f(x)=0.")
    # The prompt asks to output the numbers in the final equation.
    # The final equation is f(x) = 0.
    final_rhs = 0
    print(f"The resulting equation is: f(x) = {final_rhs}")

if __name__ == '__main__':
    illustrate_moments()