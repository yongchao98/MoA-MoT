import math

def calculate_conductance_moment(k):
    """
    Calculates the k-th moment of the dimensionless conductance 'g' for a
    disordered Majorana wire at the critical point (Class D).
    The formula is: <g^k> = Gamma(k + 0.5) / (sqrt(pi) * k!)
    """
    numerator = math.gamma(k + 0.5)
    denominator = math.sqrt(math.pi) * math.factorial(k)
    return numerator / denominator

# The goal is to find the ratio <g^4> / <g>
k_avg = 1
k_fourth = 4

# Step 1: Calculate the average conductance <g> (the first moment, k=1)
avg_g = calculate_conductance_moment(k_avg)

# Step 2: Calculate the fourth moment <g^4> (k=4)
g4_moment = calculate_conductance_moment(k_fourth)

# Step 3: Calculate the final ratio
ratio = g4_moment / avg_g

print("This script calculates the ratio between the fourth statistical moment of the dimensionless conductance and its average value.")
print("The formula for the k-th moment is: <g^k> = Gamma(k + 0.5) / (sqrt(pi) * k!)\n")

print(f"First, we calculate the average conductance <g> for k = {k_avg}:")
print(f"<g> = {avg_g}\n")

print(f"Next, we calculate the fourth moment <g^4> for k = {k_fourth}:")
print(f"<g^4> = {g4_moment}\n")

print("Finally, we compute the ratio <g^4> / <g>:")
print(f"Ratio = {g4_moment} / {avg_g} = {ratio}")
print("The exact value of this ratio is 35/64.")