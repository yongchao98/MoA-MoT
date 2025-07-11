# This script will derive the relationship between the focal length (f) and time (t)
# for a liquid-mirror telescope spun up by a constant power source.

print("Derivation of the relationship between focal length and time for a liquid mirror.")
print("-" * 75)

# Step 1: Relationship between focal length (f) and angular speed (ω)
print("Step 1: Relating focal length 'f' to angular speed 'ω'.")
print("The surface of a rotating liquid forms a paraboloid. The equation for its height 'z' at a radial distance 'x' is z = (ω^2 * x^2) / (2*g), where 'g' is the acceleration due to gravity.")
print("The equation for a standard parabolic mirror with focal length 'f' is z = x^2 / (4*f).")
print("By comparing these two equations, we find that 1 / (4*f) = ω^2 / (2*g).")
print("Solving for 'f' gives: f = g / (2 * ω^2).")
print("This shows that the focal length is inversely proportional to the square of the angular speed: f ∝ 1 / ω^2, or f ∝ ω^(-2).\n")

# Step 2: Relationship between angular speed (ω) and time (t)
print("Step 2: Relating angular speed 'ω' to time 't' for constant power 'P'.")
print("The power 'P' delivered to a rotating system is the product of torque (τ) and angular speed (ω): P = τ * ω.")
print("Torque is also given by τ = I * α, where 'I' is the moment of inertia and 'α' is the angular acceleration (dω/dt).")
print("Substituting for torque, we get the governing equation: P = I * (dω/dt) * ω.")
print("Since power 'P' and moment of inertia 'I' are constant, we can solve this differential equation by separating variables: ω dω = (P/I) dt.")
print("Integrating both sides from the initial state (ω=0 at t=0) gives: ∫ω dω = ∫(P/I) dt.")
print("The result of the integration is: (1/2)ω^2 = (P/I) * t.")
print("This shows that ω^2 is directly proportional to time: ω^2 ∝ t^1.\n")

# Step 3: Combine relationships to find f(t)
print("Step 3: Combining the relationships to find how 'f' depends on 't'.")
print("From Step 1, we have: f ∝ ω^(-2).")
print("From Step 2, we have: ω^2 ∝ t^1.")
print("Substituting the expression for ω^2 from Step 2 into the relation from Step 1 gives:")
print("f ∝ (t^1)^(-1), which simplifies to f ∝ t^(-1).\n")

# Step 4: Final Answer
print("Step 4: Determining the value of 'n'.")
print("The problem asks for the value of 'n' in the relationship f ∝ t^n.")
n = -1
print(f"Our derived relationship is f ∝ t^n, where the final equation is f ∝ t^({n}).")
print(f"Therefore, the value of n is {n}.")
