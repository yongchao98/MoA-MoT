import math

# Based on the reasoning, the problem has three possible outcomes depending on the values of a, b, and c.
# Case 1: a, b, and c are all different. The max and min values are at two different vertices. Distance = 1.
# Case 2: Exactly two of a, b, c are equal. The max/min is on an edge, and the other is at a vertex. Distance = sqrt(3)/2.
# Case 3: a = b = c. All values are the same. Distance = 0.

# The problem phrasing with distinct variables a, b, c suggests we should assume the general case
# where the values are distinct.
# In this general case, the distance is the side length of the triangle.

side_length = 1
distance = side_length

print("Assuming the general case where the numbers a, b, and c at the vertices are distinct.")
print("The maximum value will be at one vertex and the minimum value at another.")
# The distance between any two vertices of an equilateral triangle with side length 1 is 1.
print("The final equation for the distance is simply the side length.")
print(f"Distance = {distance}")

final_answer = distance
# The final answer is the number from the equation above.
print(final_answer)