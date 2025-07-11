import numpy as np
from scipy.integrate import quad

def borwein_integrals_analysis():
    """
    Calculates the Borwein integrals I_n for n=1 to 9 and prints the results
    to help analyze the properties of the sequence.
    """

    # Define a numerically stable sinc function for sin(x)/x
    def custom_sinc(x):
        # numpy.sinc(x) computes sin(pi*x)/(pi*x). We want sin(x)/x.
        # So we can use np.sinc(x/np.pi).
        return np.sinc(x / np.pi)

    # Define the integrand for I_n, which is the product of sinc functions
    def integrand(x, n):
        """
        The integrand function for the Borwein integral I_n.
        Product of sinc(x/k) from k=1 to n.
        """
        if n < 1:
            return 0.0
        
        # The product starts as 1.0
        product = 1.0
        for k in range(1, n + 1):
            product *= custom_sinc(x / k)
        return product

    # The reference value for comparison
    pi_half = np.pi / 2
    print(f"Analyzing Borwein Integrals I_n = integral from 0 to inf of product[sin(x/k)/(x/k)] dx")
    print(f"Reference value pi/2 = {pi_half:.15f}\n")

    # Loop through n from 1 to 9 to calculate and display each integral
    for n in range(1, 10):
        # Use scipy.integrate.quad for numerical integration
        # The args parameter passes 'n' to the integrand function
        result, error = quad(integrand, 0, np.inf, args=(n,))
        
        # The problem asks to output each number in the final equation.
        # We will format the output as "I_n = [value]".
        print(f"For n = {n}:")
        print(f"  I_{n} = {result:.15f}")
        print(f"  Numerical Error Estimate = {error:.2e}")
        print(f"  Difference (I_{n} - pi/2) = {result - pi_half:+.15f}")
        print("-" * 40)

# Execute the analysis
borwein_integrals_analysis()