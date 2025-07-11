import cmath

def calculate_bare_greens_function():
    """
    Calculates the bare Green's function for a given set of parameters
    and prints the steps.
    """
    # --- Parameters ---
    # Let's choose some example values for the variables.
    # In many theoretical calculations, h-bar (Planck's constant) is set to 1.
    
    # omega: The frequency of the particle/excitation
    omega = 2.5
    
    # epsilon_k: The single-particle energy eigenvalue for state k
    epsilon_k = 2.0
    
    # eta: A small positive infinitesimal to handle the pole correctly
    eta = 0.1

    # --- Calculation ---
    # The formula for the bare Green's function is G_0 = 1 / (omega - epsilon_k + i*eta)
    # In Python, the imaginary unit 'i' is represented as '1j'.
    denominator = (omega - epsilon_k) + 1j * eta
    G_0 = 1 / denominator

    # --- Output ---
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy epsilon_k is:")
    print("G_0(k, omega) = 1 / (omega - epsilon_k + i*eta)\n")
    print("This shows that G_0 is inversely proportional to the term (omega - epsilon_k).\n")
    print("-" * 50)
    print("Here is a numerical example:\n")

    print(f"Given the parameters:")
    print(f"  Frequency (omega)            = {omega}")
    print(f"  Single-particle energy (epsilon_k) = {epsilon_k}")
    print(f"  Infinitesimal (eta)          = {eta}\n")

    print("The final equation is constructed as follows:")
    # This fulfills the requirement to "output each number in the final equation"
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")
    print(f"G_0 = 1 / ({denominator})")
    print(f"G_0 = {G_0}")

    # For clarity, also print the real and imaginary parts separately
    print(f"\nThe resulting complex number is:")
    print(f"  Real Part      = {G_0.real:.4f}")
    print(f"  Imaginary Part = {G_0.imag:.4f}")

# Execute the function
calculate_bare_greens_function()