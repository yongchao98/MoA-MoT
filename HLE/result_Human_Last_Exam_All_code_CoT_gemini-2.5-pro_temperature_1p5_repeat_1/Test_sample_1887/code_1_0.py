class Ordinal:
    """A simple class to represent and print complex ordinals."""

    def __init__(self, op, args=None):
        """
        Initializes an Ordinal object.
        'op' can be 'const' for a number, 'omega' for ω_..., or 'sum' for addition.
        """
        self.op = op
        self.args = args

    def __str__(self):
        """Creates a string representation of the ordinal."""
        if self.op == 'const':
            return str(self.args)
        if self.op == 'omega':
            if self.args is None:
                return "ω"
            # Use curly braces for subscripts for clarity
            return f"ω_{{{str(self.args)}}}"
        if self.op == 'sum':
            return " + ".join(map(str, self.args))
        return ""

def solve():
    """
    Constructs and prints the ordinal ω_{ω_{ω+5}}.
    """
    # The term ω+5
    omega = Ordinal('omega', None)
    five = Ordinal('const', 5)
    omega_plus_five = Ordinal('sum', [omega, five])

    # The term ω_{ω+5}
    omega_sub_omega_plus_five = Ordinal('omega', omega_plus_five)

    # The final ordinal ω_{ω_{ω+5}}
    final_ordinal = Ordinal('omega', omega_sub_omega_plus_five)

    print("The final equation for the order type is:")
    final_equation = f"Order Type = {str(final_ordinal)}"
    print(final_equation)
    print("\nThe numbers present in this ordinal expression are:")
    print(5)

solve()