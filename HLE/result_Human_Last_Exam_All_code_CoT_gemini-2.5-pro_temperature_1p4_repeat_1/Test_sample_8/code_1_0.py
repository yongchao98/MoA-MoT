import sys
from io import StringIO

# A simple class to represent a symbolic expression like '2+alpha'
class SymbolicExpression:
    def __init__(self, numeric_part, symbolic_part):
        self.numeric = numeric_part
        self.symbolic = symbolic_part

    def __str__(self):
        if self.numeric == 0:
            return self.symbolic
        if not self.symbolic:
            return str(self.numeric)
        return f"{self.numeric} + {self.symbolic}"

def solve_conormal_space_problem():
    """
    Calculates the order of the conormal space for R(sigma)f.

    This problem involves determining the regularity of a function after applying
    the resolvent of the wave operator. The key principle is that the resolvent
    is a smoothing operator.
    """

    # 1. The initial function f is in A^(2+alpha)(X).
    #    The order 's' of this conormal space is 2+alpha.
    s = SymbolicExpression(2, "alpha")
    print(f"Step 1: The input function f belongs to the conormal space A^s(X).")
    print(f"         The order of this space is s = {s}.\n")

    # 2. The operator is the resolvent R(sigma) = (Box_g - sigma^2)^-1.
    #    The wave operator Box_g is a second-order differential operator (order 2).
    #    The resolvent, being its inverse, is a pseudodifferential operator of order m = -2.
    m = -2
    print(f"Step 2: The operator is the resolvent R(sigma).")
    print(f"         The order of this operator is m = {m}.\n")

    # 3. Apply the mapping theorem for pseudodifferential operators.
    #    An operator of order 'm' maps A^s to A^(s-m).
    #    The new order s_new = s - m.
    print(f"Step 3: According to the mapping property of pseudodifferential operators,")
    print(f"         an operator of order m maps a space of order s to a space of order s - m.\n")


    # 4. Calculate the new order.
    #    s_new = (2 + alpha) - (-2) = 4 + alpha
    s_new_numeric = s.numeric - m
    s_new = SymbolicExpression(s_new_numeric, s.symbolic)
    
    print(f"Step 4: Calculating the new order s_new:")
    print(f"         s_new = s - m")
    print(f"         s_new = ({s}) - ({m})")
    print(f"         The equation for the final order is: {s.numeric} - ({m}) = {s_new_numeric}")
    print(f"         So, the new symbolic order is s_new = {s_new}.\n")

    # 5. State the final result.
    final_space = f"A^({s_new})(X)"
    print(f"Conclusion: The function R(sigma)f belongs to the conormal space {final_space}.")

    # Storing the final answer in the required format
    # This part will not be printed to the console but captured for the final answer.
    # We use string manipulation to prepare the answer string.
    old_stdout = sys.stdout
    sys.stdout = captured_output = StringIO()
    print(f"<<<A^({s_new})(X)>>>")
    sys.stdout = old_stdout

    return captured_output.getvalue().strip()


# Run the solver and print its output
solve_conormal_space_problem()