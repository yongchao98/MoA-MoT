import sympy
from sympy import Function, Symbol, Sum, sin, cos, Eq, oo

def explain_spectral_series_for_toroid():
    """
    Explains and demonstrates the appropriate spectral series for
    poloidal dependence in toroidal systems.
    """
    # 1. The poloidal angle, here represented by 'θ' (theta), is periodic.
    #    This means any function f(θ) repeats every 2π.
    # 2. A Fourier series is the mathematical tool for expanding periodic functions
    #    into a sum of sines and cosines.
    
    # Define symbolic variables for the equation
    theta = Symbol('θ')  # The periodic poloidal angle
    n = Symbol('n', integer=True, positive=True) # The summation index (mode number)
    a_0 = Symbol('a_0')  # The constant term (average value)
    a_n = Symbol(f'a_n') # The cosine coefficients
    b_n = Symbol(f'b_n') # The sine coefficients
    f = Function('f')

    # Construct the symbolic Fourier series expansion
    # f(θ) = a_0 + Σ [ a_n*cos(n*θ) + b_n*sin(n*θ) ] for n=1 to infinity
    fourier_sum = Sum(a_n * cos(n * theta) + b_n * sin(n * theta), (n, 1, oo))
    fourier_series_equation = Eq(f(theta), a_0 + fourier_sum)

    # Print the explanation and the equation
    print("The poloidal dependence in a toroidal system is described by a periodic angle (θ).")
    print("The Fourier series is the standard spectral expansion for any periodic function.")
    print("The code below prints the general form of this series expansion:")
    print("-" * 60)
    # The pretty_print function will output the equation with all its symbols and indices ("numbers")
    sympy.pretty_print(fourier_series_equation)
    print("-" * 60)
    print("This shows that the function f(θ) is represented by a sum of trigonometric functions (sines and cosines), which is the definition of a Fourier series.")
    print("\nTherefore, Fourier series is the correct technique.")

if __name__ == '__main__':
    explain_spectral_series_for_toroid()