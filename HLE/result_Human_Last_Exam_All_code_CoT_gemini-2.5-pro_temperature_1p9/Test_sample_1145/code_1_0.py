import numpy as np
from scipy.integrate import quad

def main():
    """
    Calculates the probability of finding a particle in a 1D box in a given interval.
    """
    # 1. Define the constants from the problem
    n = 2           # Quantum state
    u1 = 0.495      # Lower bound of the interval (as a fraction of box length a)
    u2 = 0.505      # Upper bound of the interval (as a fraction of box length a)

    # 2. Define the integrand function for the probability
    # This corresponds to 2 * sin^2(n * pi * u)
    def probability_density(u):
        return 2 * (np.sin(n * np.pi * u))**2

    # 3. Perform the numerical integration using scipy.integrate.quad
    # The quad function returns the result of the integral and an estimate of the error
    probability, error = quad(probability_density, u1, u2)
    
    # 4. Print the final results, showing the components of the calculation as requested
    print("--- Calculation of Probability for a Particle in a 1D Box ---")
    print(f"The particle is in the quantum state n = {n}.")
    print(f"We are finding the probability in the interval from x = {u1}*a to x = {u2}*a.")
    print("The equation for the probability is the integral of P(u) = 2 * sin^2(n*pi*u) from u1 to u2.")
    print(f"P = Integral from {u1} to {u2} of 2 * sin^2({n}*pi*u) du")
    print(f"\nThe calculated probability is: {probability}")
    # Note: The interval is centered on x=0.5a, which is a node for the n=2 state
    # (where the wavefunction is zero). Therefore, the probability is expected to be very small.

if __name__ == "__main__":
    main()