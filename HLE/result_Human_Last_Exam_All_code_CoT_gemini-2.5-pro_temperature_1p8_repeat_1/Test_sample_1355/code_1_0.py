import math

def calculate_conductance_moment(n: int) -> float:
    """
    Calculates the n-th statistical moment of the dimensionless conductance 'g'
    for a disordered Majorana wire at the critical point.
    The formula is <g^n> = Γ(n + 1/2) / (sqrt(π) * n!).
    """
    if n < 0 or not isinstance(n, int):
        raise ValueError("Moment 'n' must be a non-negative integer.")
    
    # Numerator: Gamma function Γ(n + 1/2)
    numerator = math.gamma(n + 0.5)
    
    # Denominator: sqrt(π) * n!
    denominator = math.sqrt(math.pi) * math.factorial(n)
    
    return numerator / denominator

# The user wants to find the ratio of the fourth moment to the average value.
# The average value is the first moment (n=1).
n_avg = 1
# The fourth moment is n=4.
n_fourth = 4

# Calculate the average value <g>
avg_g = calculate_conductance_moment(n_avg)

# Calculate the fourth moment <g^4>
fourth_moment_g = calculate_conductance_moment(n_fourth)

# Calculate the required ratio
if avg_g == 0:
    ratio = float('inf')
else:
    ratio = fourth_moment_g / avg_g

# As per the request, we output each number in the final equation.
print("This problem concerns the universal conductance statistics at the critical point of disordered Majorana wires.")
print("The average dimensionless conductance <g> is the 1st moment (n=1).")
print(f"Value of <g>: {avg_g}")
print("\nThe fourth statistical moment <g^4> corresponds to n=4.")
print(f"Value of <g^4>: {fourth_moment_g}")
print("\nThe required ratio is <g^4> / <g>.")
print(f"Final Equation: {fourth_moment_g} / {avg_g} = {ratio}")
print("\nIn fractional form, <g^4> = 35/128, <g> = 1/2, so the ratio is (35/128) / (1/2) = 35/64.")
