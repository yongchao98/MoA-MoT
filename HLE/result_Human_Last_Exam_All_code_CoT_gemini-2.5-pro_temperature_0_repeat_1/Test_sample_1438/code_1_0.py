# A simple class to represent a Grassmann expression of the form: c0 + c1*η
class GrassmannExpr:
    """
    Represents a simple polynomial in a single Grassmann variable η.
    Since η^2 = 0, any expression can be written as c0 + c1*η.
    """
    def __init__(self, c0, c1):
        # c0 is the constant part, c1 is the coefficient of η
        self.c0 = c0
        self.c1 = c1

    def __repr__(self):
        # A user-friendly string representation of the expression
        if self.c0 == 0 and self.c1 == 0:
            return "0"
        
        parts = []
        if self.c0 != 0:
            parts.append(str(self.c0))
        
        if self.c1 != 0:
            if self.c1 == 1:
                # Avoid printing 1η, just print η
                parts.append("η")
            elif self.c1 == -1:
                # Avoid printing + -1η, just print - η
                parts.append("-η")
            else:
                # Format the coefficient with the variable
                parts.append(f"{self.c1}η")
        
        # Join the parts with a plus sign and clean up formatting
        return " + ".join(parts).replace(" + -", " - ")

def berezin_integral(expr):
    """
    Applies the Berezin integration rules for a single variable:
    ∫ dη = 0
    ∫ η dη = 1
    This means the integral of an expression (c0 + c1*η) is simply c1.
    """
    # The integral "picks out" the coefficient of η
    return expr.c1

# --- Main Demonstration ---

# Represent the Grassmann variable η itself.
# In our class notation, η is equivalent to the expression 0 + 1*η.
eta_variable = GrassmannExpr(c0=0, c1=1)

# Perform the integral ∫ η dη using our defined function.
# The function berezin_integral(eta_variable) will return the coefficient of η, which is 1.
result = berezin_integral(eta_variable)

# Print the final equation as requested.
# The integral symbol ∫ is represented by 'Integral'.
# The measure dη is implied by the integration function.
print("The defining normalization for the Grassmann integral measure is:")
print(f"Integral( {eta_variable} dη ) = {result}")
