import math

def solve_expected_value():
    """
    Calculates and explains the expected value of T.
    
    The problem is to find E[T] where T = sum_{i=1 to m} Y_i / 3^i.
    Y_i are i.i.d U(0,1) and m is the smallest integer such that Y_m > Y_{m+1}.

    Step 1: Express E[T] as an infinite sum.
    E[T] = E[sum_{i=1 to inf} (Y_i / 3^i) * I(i <= m)]
         = sum_{i=1 to inf} (1 / 3^i) * E[Y_i * I(i <= m)]
    
    Step 2: Evaluate the expectation term.
    The event {i <= m} is equivalent to {Y_1 <= Y_2 <= ... <= Y_i}. Let's call this event B_i.
    E[Y_i * I(B_i)] = P(B_i) * E[Y_i | B_i].
    
    P(B_i) is the probability of one specific ordering of i i.i.d. variables, so P(B_i) = 1/i!.
    E[Y_i | B_i] is the expected value of the maximum of i U(0,1) variables, which is i/(i+1).
    So, E[Y_i * I(B_i)] = (1/i!) * (i/(i+1)) = i / ((i+1)!).

    Step 3: Sum the series.
    E[T] = sum_{i=1 to inf} (1 / 3^i) * (i / ((i+1)!))
    We use the identity i/((i+1)!) = 1/i! - 1/((i+1)!).
    E[T] = sum_{i=1 to inf} (1/3^i) * (1/i! - 1/((i+1)!))
         = sum_{i=1 to inf} (1/3)^i / i! - sum_{i=1 to inf} (1/3)^i / ((i+1)!)
    
    The first sum is (e^(1/3) - 1).
    The second sum is 3 * (e^(1/3) - 1 - 1/3) = 3*e^(1/3) - 4.
    
    Step 4: Final result.
    E[T] = (e^(1/3) - 1) - (3*e^(1/3) - 4)
         = e^(1/3) - 1 - 3*e^(1/3) + 4
         = 3 - 2*e^(1/3)
    """
    
    # Final equation is E[T] = 3 - 2 * e^(1/3)
    c1 = 3
    c2 = 2
    exponent = 1/3
    
    # Calculate the numerical value
    exp_val = math.exp(exponent)
    result = c1 - c2 * exp_val
    
    print("The derivation leads to the final expression for the expected value E[T]:")
    print("\nE[T] = 3 - 2 * e^(1/3)\n")
    print("Here is the step-by-step calculation of the numerical value:")
    print(f"Let c1 = {c1}")
    print(f"Let c2 = {c2}")
    print(f"Let x = 1/3 = {exponent:.6f}")
    print(f"Then e^x = exp({exponent:.6f}) = {exp_val:.6f}")
    print(f"E[T] = c1 - c2 * e^x = {c1} - {c2} * {exp_val:.6f} = {result:.6f}")
    
    return result

# Execute the function to print the solution.
solve_expected_value()
