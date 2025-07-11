class GrassmannFunction:
    """
    A simple class to represent a function of a single Grassmann variable psi.
    Since psi*psi = 0, any function can be written as f(psi) = A + B*psi.
    We represent this as an object with attributes self.const = A and self.linear = B.
    """
    def __init__(self, const_part, linear_part):
        self.const = const_part  # The constant part (A)
        self.linear = linear_part # The part proportional to psi (B)

    def __str__(self):
        # Creates a string representation like "A + B*psi"
        const_str = str(self.const)
        linear_str = str(self.linear)

        if self.const == 0 and self.linear == 0:
            return "0"
        if self.linear == 0:
            return const_str
        if self.const == 0:
            if self.linear == 1:
                return "psi"
            elif self.linear == -1:
                return "-psi"
            else:
                return f"{linear_str}*psi"
        
        op = "+" if self.linear > 0 else "-"
        if abs(self.linear) == 1:
             return f"{const_str} {op} psi"
        else:
             return f"{const_str} {op} {abs(self.linear)}*psi"
             
    def __mul__(self, other):
        # Implements (A + B*psi) * (C + D*psi) = A*C + (A*D + B*C)*psi
        # The B*D*psi*psi term is zero.
        new_const = self.const * other.const
        new_linear = self.const * other.linear + self.linear * other.const
        return GrassmannFunction(new_const, new_linear)

def berezin_integral(func):
    """
    Applies the Berezin integration rules:
    ∫ d(psi) (A + B*psi) = B
    The integral of a Grassmann function is defined to be the coefficient
    of the linear term.
    """
    return func.linear

# --- Main Demonstration ---

# Define the base Grassmann variable, psi = 0 + 1*psi
psi = GrassmannFunction(0, 1)

# Define a constant function, one = 1 + 0*psi
one = GrassmannFunction(1, 0)

print("In the Grassmann algebra, the Pauli exclusion principle is encoded by the property psi^2 = 0.")
print(f"Let's represent our variable as: psi = {psi}")
psi_squared = psi * psi
print(f"Calculating psi * psi results in: {psi_squared}")
print("-" * 30)

print("The integral measure 'd(psi)' is defined by the following rules:")
print("Rule 1: The integral of a constant is zero.")
print("Rule 2: The integral of the variable itself is one.")
print("\nLet's demonstrate these rules.")

# 1. Integrate the constant function '1'
integral_of_one = berezin_integral(one)
print(f"\nThe function to integrate is f(psi) = {one}")
print(f"The result of the integral ∫ d(psi) * ({one}) is: {integral_of_one}")
print(f"So, the final equation is: ∫ d(psi) * 1 = {integral_of_one}")


# 2. Integrate the function 'psi'
integral_of_psi = berezin_integral(psi)
print(f"\nThe function to integrate is f(psi) = {psi}")
print(f"The result of the integral ∫ d(psi) * ({psi}) is: {integral_of_psi}")
print(f"So, the final equation is: ∫ d(psi) * psi = {integral_of_psi}")
