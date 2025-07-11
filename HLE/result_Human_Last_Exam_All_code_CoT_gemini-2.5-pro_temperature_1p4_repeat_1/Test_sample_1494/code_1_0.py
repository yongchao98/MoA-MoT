import math

def calculate_decay(t, C=100.0, power=-1.5):
  """
  Calculates the value of a quantity decaying algebraically with time.
  The decay rate is t^power.
  
  Args:
    t: The time.
    C: A constant prefactor.
    power: The exponent of time t.
  
  Returns:
    The decayed value.
  """
  return C * (t**power)

# We are modeling the decay of ||∂_x ρ(t)||_L2, which is proportional to t^(-3/2).
# Let's represent this relationship with an equation y = C * t^(-3/2).
# Here, y represents ||∂_x ρ(t)||_L2 and -3/2 is the power of decay.
C = 100.0
power = -1.5 # This corresponds to -3/2

# Time points for calculation
times = [1, 2, 5, 10]

# The decay equation is y = C * t^(power)
print(f"The decay is modeled by the equation: y = {C} * t^({power})")
print("This corresponds to ||∂_x ρ(t)||_L2 ∝ 1 / t^(3/2)")
print("\nLet's calculate the values at different times t:")

# Calculate and print the result for each time point, showing the full equation.
for t in times:
    result = calculate_decay(t, C, power)
    # Print the equation with numbers plugged in.
    # We explicitly show the power as 3/2 in the denominator for clarity.
    denominator_val = t ** 1.5
    print(f"For t = {t}: y = {C} * (1 / {t}^(3/2)) = {C} / {denominator_val:.3f} = {result:.3f}")
