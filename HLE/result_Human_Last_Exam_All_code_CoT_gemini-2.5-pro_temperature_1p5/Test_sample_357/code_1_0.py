import math

# For the Low type contract (v_L, p_L)
v_L = 2/3
p_L = (v_L)**2 - (1/3)*v_L + 20

# For the High type contract (v_H, p_H)
# We solve the quadratic equation v_H^2 - v_H - 92/9 = 0 for v_H
v_H = (1 + math.sqrt(1 + 4 * (92/9))) / 2
p_H = v_H + 30

# Print the contracts in the final format
# We format them into an equation as requested
# Each number is outputted
print(f"The separating equilibrium contracts are:")
print(f"(v_L, p_L) = ({v_L:.3f}, {p_L:.3f})")
print(f"(v_H, p_H) = ({v_H:.3f}, {p_H:.3f})")

# Let's write the final answer as a single string equation
final_equation = f"(v_L, p_L) = ({v_L}, {p_L}), (v_H, p_H) = ({v_H}, {p_H})"
# Although the prompt asks for print, this is the final answer string I will wrap
# For better readability, I will also print the exact fractional forms where simple.
# v_L = 2/3, p_L = 182/9
v_L_f = "2/3"
p_L_f = "182/9"
# v_H and p_H are irrational
v_H_f = f"(3 + sqrt(377))/6"
p_H_f = f"30 + (3 + sqrt(377))/6"

print("\nWith exact values:")
print(f"(v_L, p_L) = ({v_L_f}, {p_L_f})")
print(f"(v_H, p_H) = ({v_H_f}, {p_H_f})")
