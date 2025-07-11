import math

def solve_expected_value():
    """
    This function explains the derivation and calculates the expected value of T.
    """
    print("Step 1: Express the expected value E[T] as an infinite sum.")
    print("E[T] = E[Sum_{i=1 to m} Y_i/3^i]")
    print("Using the identity for expectations of stopped sums with non-negative terms, we get:")
    print("E[T] = Sum_{i=1 to inf} E[ (Y_i/3^i) * I(m >= i) ]")
    print("where I(m >= i) is the indicator function for the event m >= i.")
    print("-" * 20)

    print("Step 2: Determine the condition for m >= i.")
    print("The event {m >= i} means Y_1 <= Y_2 <= ... <= Y_i.")
    print("So, E[T] = Sum_{i=1 to inf} (1/3^i) * E[Y_i * I(Y_1 <= ... <= Y_i)]")
    print("-" * 20)

    print("Step 3: Calculate the expectation term E[Y_i * I(Y_1 <= ... <= Y_i)].")
    print("This is given by the integral of y_i over the region 0 <= y_1 <= ... <= y_i <= 1.")
    print("Integral = integral from 0 to 1 of y_i * (y_i^(i-1) / (i-1)!) dy_i")
    print("= 1/((i-1)!) * integral from 0 to 1 of y_i^i dy_i")
    print("= 1/((i-1)!) * (1/(i+1)) = i / (i+1)!")
    print("-" * 20)

    print("Step 4: Substitute this back into the sum for E[T].")
    print("E[T] = Sum_{i=1 to inf} (1/3^i) * (i / (i+1)!)")
    print("-" * 20)

    print("Step 5: Evaluate the infinite series.")
    print("We can rewrite the term i / (i+1)! as (i+1 - 1) / (i+1)! = 1/i! - 1/(i+1)!")
    print("So, E[T] = Sum_{i=1 to inf} [ 1/(3^i * i!) - 1/(3^i * (i+1)!) ]")
    print("This can be split into two sums:")
    print("Sum 1: Sum_{i=1 to inf} (1/3)^i / i!")
    print("Sum 2: Sum_{i=1 to inf} 1 / (3^i * (i+1)!)")
    print("\nEvaluating Sum 1:")
    print("The Taylor series for e^x is Sum_{i=0 to inf} x^i / i!.")
    print("Sum 1 is e^(1/3) - (1/3)^0/0! = e^(1/3) - 1.")
    sum1_val = math.exp(1/3) - 1
    
    print("\nEvaluating Sum 2:")
    print("Let j = i+1. The sum becomes Sum_{j=2 to inf} 1 / (3^(j-1) * j!)")
    print("= 3 * Sum_{j=2 to inf} (1/3)^j / j!")
    print("= 3 * (e^(1/3) - (1/3)^0/0! - (1/3)^1/1!)")
    print("= 3 * (e^(1/3) - 1 - 1/3) = 3 * (e^(1/3) - 4/3) = 3*e^(1/3) - 4.")
    sum2_val = 3 * math.exp(1/3) - 4
    print("-" * 20)
    
    print("Step 6: Combine the results to find the final answer.")
    print("E[T] = (Sum 1) - (Sum 2)")
    print("E[T] = (e^(1/3) - 1) - (3*e^(1/3) - 4)")
    print("E[T] = e^(1/3) - 1 - 3*e^(1/3) + 4")
    final_symbolic_expr = "3 - 2*e^(1/3)"
    print(f"E[T] = {final_symbolic_expr}")
    
    # As requested, output each number in the final equation.
    print("\nThe final equation is E[T] = a - b * exp(c/d), where:")
    a = 3
    b = 2
    c = 1
    d = 3
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"d = {d}")
    print(f"So, E[T] = {a} - {b} * e^({c}/{d})")

    # Calculate the numerical value
    final_numerical_value = a - b * math.exp(c/d)
    print(f"\nThe numerical value is approximately: {final_numerical_value}")

if __name__ == '__main__':
    solve_expected_value()