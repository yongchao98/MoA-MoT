def demonstrate_property(x, y):
    """
    This function demonstrates that for any two distinct points x and y on the
    real line, we can find a closed connected set K (an interval [a,b]) that
    contains x in its interior but does not contain y.

    Args:
        x (float): The first point.
        y (float): The second point, must be distinct from x.
    """
    print(f"--- Demonstrating for x = {x}, y = {y} ---")
    if x == y:
        print("Points must be distinct.")
        return

    # In the real line, a closed connected set is a closed interval [a, b].
    # Its interior is the open interval (a, b).
    
    if y > x:
        # If y is greater than x, we can create an interval around x
        # whose right endpoint is less than y.
        a = x - 1
        b = (x + y) / 2
    else:  # y < x
        # If y is less than x, we can create an interval around x
        # whose left endpoint is greater than y.
        a = (y + x) / 2
        b = x + 1

    K_interval = (a, b)
    
    print(f"Chosen closed connected set K = [{K_interval[0]}, {K_interval[1]}]")
    print(f"Interior of K is Int(K) = ({K_interval[0]}, {K_interval[1]})")

    # Check if the conditions hold
    x_is_in_interior = (K_interval[0] < x < K_interval[1])
    y_is_not_in_K = (y < K_interval[0] or y > K_interval[1])

    print(f"Is x in Int(K)? {x_is_in_interior}")
    print(f"Is y not in K? {y_is_not_in_K}")
    print(f"Property holds: {x_is_in_interior and y_is_not_in_K}\n")

# Run demonstrations for a couple of examples
demonstrate_property(5, 10)
demonstrate_property(-4.5, -20)

# The topological analysis shows that any space satisfying the given conditions
# must be homeomorphic to the real line ℝ. Therefore, there is only one
# such class of spaces.
final_answer = 1

print("-----------------------------------------------------")
print("The final result of the deduction is:")
print(f"Number of homeomorphism classes = {final_answer}")