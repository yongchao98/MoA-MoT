import numpy as np

def illustrate_wasserstein_gradient():
    """
    Illustrates that the Wasserstein gradient of J(mu) = 0.5 * W(mu, nu)^2
    is the trivial vector at the minimum mu = nu.
    """

    print("--- Theoretical Setup for 1D Gaussian Case ---")
    print("Let the target measure nu be the standard normal distribution N(mean=0, std=1).")
    print("Let mu_a be a test measure from the family N(mean=a, std=1).")
    print("The functional is J(mu_a) = 0.5 * W(mu_a, nu)^2.")
    print("The minimum is at mu = nu, which corresponds to a = 0.")
    print("\nIn this Gaussian case, we have simple analytical formulas:")
    # Formula for W^2 between N(m1, s1^2) and N(m2, s2^2) is (m1-m2)^2 + (s1-s2)^2
    # Here, m1=a, m2=0, s1=1, s2=1.
    print("  1. Wasserstein distance squared: W(mu_a, nu)^2 = (a - 0)^2 + (1 - 1)^2 = a^2")
    print("  2. The functional J: J(mu_a) = 0.5 * a^2")
    # Optimal map from N(a, 1) to N(0, 1) is T(x) = x - a.
    # The gradient vector field is v(x) = x - T(x).
    print("  3. The gradient vector field v(x) = x - (x - a) = a")
    print("--- End of Setup ---\n")

    def compute_and_show_gradient_at_a(a):
        """
        Computes and displays the gradient for a given parameter 'a'.
        """
        # The gradient vector field is the constant value 'a'.
        gradient_value = a
        # The value of the functional J is 0.5 * a^2
        j_value = 0.5 * a**2

        print(f"--- Analysis for a = {a} ---")
        print(f"The measure is mu_{a} = N(mean={a}, std=1)")
        # Final equation output for J
        print(f"The functional value is J(mu_{a}) = 0.5 * {a}^2 = {j_value}")
        # Final equation output for the gradient
        print(f"The gradient vector field v(x) is the constant function v(x) = {gradient_value}")
        
        if a == 0:
            print("\nThis is the minimum of the functional J (since mu_0 = nu).")
            print("The gradient vector field v(x) = 0 is the 'trivial tangent vector'.")
        else:
            print("\nThis is NOT the minimum of the functional J.")
            print("The gradient vector field is non-trivial.")

    # Let's test points away from, near, and at the minimum.
    print(">>> Testing a point AWAY from the minimum:")
    compute_and_show_gradient_at_a(a=2.0)
    print("-" * 40)
    
    print(">>> Testing a point NEAR the minimum:")
    compute_and_show_gradient_at_a(a=0.01)
    print("-" * 40)
    
    print(">>> Testing AT the minimum:")
    compute_and_show_gradient_at_a(a=0.0)
    print("-" * 40)

    print("\nAs the code illustrates, the gradient systematically approaches 0 as 'a' approaches 0,")
    print("and is exactly the zero vector field when 'a' is 0, which is at the minimum of J.")

# Execute the illustration
illustrate_wasserstein_gradient()