import numpy as np
from scipy.integrate import quad

def demonstrate_scaling():
    """
    This function numerically demonstrates that the integral of a periodic
    function over a 3D ball scales with R^3.
    """
    # Amplitude of the periodic function. Let's use M=0.5.
    M = 0.5
    
    # The periodic function is f(x) = (M*cos(x))^2.
    def f(x):
        return (M * np.cos(x))**2

    # The average value of f(x) over one period is C = M^2 / 2.
    C = M**2 / 2
    
    # The integrand for the 1D integral is f(x) * (area of disk at x).
    def integrand(x, R):
        return f(x) * np.pi * (R**2 - x**2)

    # We expect the integral I(R) to be asymptotically C * (4/3)*pi*R^3.
    expected_constant = C * 4/3 * np.pi
    
    print("Demonstrating the scaling of the integral I(R) = integral_{B_R} (v'(x))^2 dV.")
    print(f"We use a proxy function |v'|^2 = ({M}*cos(x))^2.")
    print(f"The theoretical average value of |v'|^2 is C = {C:.4f}.")
    print(f"The integral I(R) should scale like K * R^a.")
    print(f"Our analysis predicts a = 3 and K = 4/3 * pi * C.")
    print(f"Let's check the ratio I(R) / R^3 for large R.\n")

    print(f"Theoretical constant K = 4/3 * pi * {C:.4f} = {expected_constant:.4f}")
    print("-" * 50)
    print("{:<10} {:<20} {:<20}".format("R", "I(R)", "I(R) / R^3"))
    print("-" * 50)
    
    for R_val in [10, 50, 100, 200, 500]:
        # Perform numerical integration from -R to R.
        integral_val, _ = quad(integrand, -R_val, R_val, args=(R_val,), limit=200)
        ratio = integral_val / (R_val**3)
        print("{:<10.1f} {:<20.4e} {:<20.4f}".format(R_val, integral_val, ratio))
    
    print("-" * 50)
    print("As R increases, the ratio I(R) / R^3 converges to the predicted constant.")
    print("\nThis supports the conclusion that the integral grows as R^3, making a=3.")
    # Outputting the numbers from the "final equation" I(R) ~ K * R^a
    print("\nFinal equation analysis:")
    print(f"The power 'a' is: 3")
    print(f"The constant 'K' for M={M} is approx: {expected_constant:.4f}")


demonstrate_scaling()