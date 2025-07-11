import math

def solve_mixing_problem():
    """
    Calculates the number of repetitions needed to make the liquids in two cups
    "the same" based on the given epsilon and delta parameters.
    """
    # Fix parameters epsilon and delta from the problem description.
    # You can change these values to test other scenarios.
    epsilon = 0.01
    delta = 0.1

    print(f"Given parameters: ε = {epsilon}, δ = {delta}")
    print("-" * 40)

    # The condition to meet is |(1 - 2*δ)^n| <= ε.
    # We solve for n: n >= ln(ε) / ln|1 - 2*δ|.
    # As n must be an integer, we need the ceiling of this value.
    
    # Handle the special case where δ = 0.5
    if delta == 0.5:
        # If delta is 0.5, 1 - 2*delta = 0. The mixture becomes perfect in one step.
        # The difference becomes 0, which is less than any positive epsilon.
        n = 1
        print("Special case: δ = 0.5. The liquids mix perfectly in 1 step.")
    
    else:
        # The final equation is n = ceil(ln(ε) / ln|1 - 2*δ|)
        print("Final equation: n = ceil(ln(ε) / ln|1 - 2*δ|)")
        print("\nCalculating the components of the equation:")

        # Calculate each part of the formula to show the steps
        numerator = math.log(epsilon)
        base = abs(1 - 2 * delta)
        denominator = math.log(base)
        
        print(f"1. Numerator: ln(ε) = ln({epsilon:.4f}) = {numerator:.4f}")
        print(f"2. Base of exponent: |1 - 2*δ| = |1 - 2*{delta}| = {base:.4f}")
        print(f"3. Denominator: ln(|1 - 2*δ|) = ln({base:.4f}) = {denominator:.4f}")
        
        # The ratio
        ratio = numerator / denominator
        print(f"4. Ratio: {numerator:.4f} / {denominator:.4f} = {ratio:.4f}")
        
        # n is the ceiling of this ratio
        n = math.ceil(ratio)
    
    print("-" * 40)
    print(f"The number of repetitions required is: {int(n)}")

solve_mixing_problem()