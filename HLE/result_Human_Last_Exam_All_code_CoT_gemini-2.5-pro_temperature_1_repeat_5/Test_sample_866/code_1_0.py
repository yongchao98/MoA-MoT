import math

# The input data representing the crease pattern at a single vertex.
# Format: [angle1, crease1_type, angle2, crease2_type, ...]
data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

# 1. Parse the input into separate lists for angles and crease assignments.
angles = [item for item in data if isinstance(item, int)]
creases = [item for item in data if isinstance(item, str)]

print("Step-by-step analysis of the crease pattern:")
print(f"Angles: {angles}")
print(f"Creases: {creases}\n")

# 2. Check Kawasaki's Theorem (the angle condition).
# This is a necessary condition for any vertex to be flat-foldable.
print("Step 1: Check Kawasaki's Theorem (Angle Condition)")

odd_indices = range(0, len(angles), 2)
even_indices = range(1, len(angles), 2)

sum_odd = sum(angles[i] for i in odd_indices)
sum_even = sum(angles[i] for i in even_indices)

# Construct the equation strings to show the calculation.
odd_angles_str = " + ".join(str(angles[i]) for i in odd_indices)
even_angles_str = " + ".join(str(angles[i]) for i in even_indices)

# Print the evaluation of the theorem.
print("The theorem states that the sum of alternating angles must be equal.")
print(f"Sum of odd-indexed angles: {odd_angles_str} = {sum_odd}")
print(f"Sum of even-indexed angles: {even_angles_str} = {sum_even}")

# 3. Determine the result based on the angle check.
# If the sums are not equal, the pattern is not flat-foldable.
final_answer = 0
if sum_odd != sum_even:
    print(f"\nResult: The sums are not equal ({sum_odd} != {sum_even}).")
    print("Kawasaki's Theorem is not satisfied.")
    print("Therefore, there are no valid assignments that can make this pattern flat-foldable.")
    final_answer = 0
else:
    # This block handles the case where the angle condition is met.
    # It would then check Maekawa's theorem.
    # For the given input, this block will not be executed.
    print("\nResult: The sums are equal. Kawasaki's Theorem is satisfied.")
    print("\nStep 2: Check Maekawa's Theorem (Crease Count Condition)")
    N = len(creases)
    print("The theorem states |Number of Mountain folds - Number of Valley folds| = 2.")
    print("This implies the total number of creases must be even.")

    if N % 2 != 0:
        print(f"The number of creases is {N}, which is odd.")
        print("Maekawa's Theorem cannot be satisfied.")
        final_answer = 0
    else:
        # This block calculates the number of combinations if both conditions are met.
        num_M_known = creases.count('M')
        num_V_known = creases.count('V')
        num_Q = creases.count('?')

        def combinations(n, k):
            if k < 0 or k > n: return 0
            return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

        total_assignments = 0
        
        # Case 1: N_M - N_V = 2
        # Solves for k, the number of '?' to become 'M'.
        k1_num = 2 - num_M_known + num_V_known + num_Q
        if k1_num >= 0 and k1_num % 2 == 0:
            k1 = k1_num // 2
            if 0 <= k1 <= num_Q:
                total_assignments += combinations(num_Q, k1)
        
        # Case 2: N_V - N_M = 2
        k2_num = num_V_known - num_M_known + num_Q - 2
        if k2_num >= 0 and k2_num % 2 == 0:
            k2 = k2_num // 2
            if 0 <= k2 <= num_Q:
                total_assignments += combinations(num_Q, k2)
        
        final_answer = total_assignments

# Print the final calculated number of different assignments.
print(f"\nTotal number of different assignments: {final_answer}")