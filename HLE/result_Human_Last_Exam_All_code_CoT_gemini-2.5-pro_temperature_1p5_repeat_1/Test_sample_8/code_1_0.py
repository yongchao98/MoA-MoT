import sys
import io

# A class to emulate the symbolic variable alpha
class SymbolicAlpha:
    def __repr__(self):
        return "alpha"

# Define the symbolic variable
alpha = SymbolicAlpha()

# 1. Define the initial parameters based on the problem statement.
# The function f is in the conormal space A^{s}(X) where s = 2 + alpha.
initial_order_numeric_part = 2

# The resolvent R(sigma) corresponds to the inverse of the second-order
# wave operator Box_g. In the relevant pseudodifferential calculus,
# this makes R(sigma) an operator of order m = -2.
operator_order = -2

# 2. Explain the mathematical principle.
print("The problem is to determine the conormal space for the function R(sigma)f.")
print("This can be solved by understanding the mapping properties of the resolvent operator.")
print("-" * 60)
print("1. The input function f is in the conormal space A^{s}(X) with order:")
print(f"   s = {initial_order_numeric_part} + {alpha}")
print("")
print("2. The resolvent R(sigma) of a second-order differential operator (like the wave operator) is a pseudodifferential operator of order m.")
print(f"   The order of the resolvent operator is m = {operator_order}.")
print("")
print("3. An operator of order m maps a conormal space of order s to one of order s - m.")
print("   Therefore, the new order s_new is calculated as:")
print("")

# 3. Perform and display the calculation of the new order.
# The "equation" is s_new = s - m = (2 + alpha) - (-2)
final_order_numeric_part = initial_order_numeric_part - operator_order

# We print each number involved in the final calculation.
print(f"   s_new = s - m")
print(f"   s_new = ({initial_order_numeric_part} + {alpha}) - ({operator_order})")
print(f"   s_new = {initial_order_numeric_part} - ({operator_order}) + {alpha}")
print(f"   s_new = {final_order_numeric_part} + {alpha}")
print("-" * 60)

# 4. State the final conclusion.
print(f"Thus, the function R(sigma)f belongs to the conormal space A^{{{final_order_numeric_part}+{alpha}}}(X).")

# Final answer in the specified format
# The question asks for the conormal space, so we provide it.
final_answer = f"A^({final_order_numeric_part}+{alpha})(X)"
sys.stdout = io.StringIO()
print(f'<<<{final_answer}>>>')
sys.stdout = sys.__stdout__
print(f'<<<{final_answer}>>>')
