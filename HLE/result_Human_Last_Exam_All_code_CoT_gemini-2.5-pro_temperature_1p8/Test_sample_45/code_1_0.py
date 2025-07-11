# This script derives the exponent 'n' in the relationship f ∝ t^n
# for a liquid-mirror telescope spun up by a constant power source.

# Step 1: Relate focal length (f) to angular velocity (ω).
print("Step 1: Finding the relationship between focal length (f) and angular velocity (ω).")
print("The surface of a rotating liquid forms a paraboloid with the equation z = (ω^2 * r^2) / (2*g).")
print("The standard equation for a parabola's focus is z = r^2 / (4*f).")
print("By comparing these two equations, we can solve for f:")
print("1 / (4*f) = ω^2 / (2*g)  =>  f = g / (2*ω^2)")
print("This gives the proportionality: f ∝ 1 / ω^2")
print("-" * 60)

# Step 2: Relate angular velocity (ω) to time (t).
print("Step 2: Finding the relationship between angular velocity (ω) and time (t).")
print("A constant power source (P) provides energy, so P = dK/dt, where K is kinetic energy.")
print("The rotational kinetic energy is K = (1/2) * I * ω^2, where I is the moment of inertia.")
print("Therefore, P = d/dt[(1/2) * I * ω^2] = (I/2) * d(ω^2)/dt.")
print("Since P and I are constant, d(ω^2)/dt is also constant. Integrating with ω(0)=0 gives:")
print("ω^2 = (2*P/I) * t")
print("This gives the proportionality: ω^2 ∝ t")
print("-" * 60)

# Step 3: Combine the relationships to find how f depends on t.
print("Step 3: Combining the two relationships to find f(t).")
print("We found that f ∝ 1 / ω^2 and ω^2 ∝ t.")
print("Substituting the expression for ω^2 into the expression for f, we get:")
print("f ∝ 1 / t, which can also be written as f ∝ t^-1.")
print("-" * 60)

# Step 4: Determine the final value of n.
print("Step 4: Final determination of the exponent n.")
print("The problem requires us to find n, where f ∝ t^n.")
print("From our derivation, we found the relationship is f ∝ t^-1.")

# Extract the final number as requested.
final_exponent = -1

print("\nThe final equation is f ∝ t^n, where n is:")
print(f"n = {final_exponent}")
print("\nSo the final equation form is:")
print(f"f ∝ t ** ({final_exponent})")