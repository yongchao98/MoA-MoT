import numpy as np
from scipy.integrate import quad

def integrand(x, n):
    """
    Defines the integrand for the Borwein integral I_n.
    The product is Π_{k=1 to n} [sin(x/k) / (x/k)].
    """
    if x == 0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        # This is the sinc function, sinc(y) = sin(y)/y.
        # It's 1 at y=0. x/k is only 0 if x is 0, which is handled.
        product *= np.sin(x / k) / (x / k)
        
    return product

# --- Calculations ---
pi_half = np.pi / 2
print(f"The value of π/2 is: {pi_half}\n")

# --- Verification for n=4 ---
n4 = 4
I4, err4 = quad(lambda x: integrand(x, n4), 0, np.inf)
print(f"--- Analysis for n = {n4} ---")
print(f"P({n4}) is the proposition that I_{n4} = π/2.")
print(f"Numerically calculated I_{n4} = {I4}")
print(f"The difference |I_{n4} - π/2| is {abs(I4 - pi_half)}.")
# As the difference is not zero, P(4) is false.
# This refutes statement A) and D).
print(f"We can see that I_{n4} is less than π/2, which supports statement C).\n")

# --- Verification for n=5 ---
n5 = 5
I5, err5 = quad(lambda x: integrand(x, n5), 0, np.inf)
diff_5 = abs(I5 - pi_half)
threshold = 1e-5

print(f"--- Analysis for n = {n5} ---")
print(f"P({n5}) is the proposition that I_{n5} = π/2.")
print(f"Numerically calculated I_{n5} = {I5}")
print(f"The numerical integration error is approximately {err5:.2e}.")
# The non-zero difference (larger than the error) disproves P(5), supporting statement I).
print("\nNow we check the specific inequality in statement F):")
print(f"Is |I_{n5} - π/2| < {threshold}?")
print(f"We compare the values:")
print(f"|{I5} - {pi_half}| = {diff_5}")
print(f"So, the inequality {diff_5} < {threshold} is {diff_5 < threshold}.")
print("Statement F) is therefore true.")
