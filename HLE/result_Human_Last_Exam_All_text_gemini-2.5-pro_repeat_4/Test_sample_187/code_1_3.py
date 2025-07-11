26     def is_right_triangle(self):
27         # Calculate squared lengths of all sides
28         sides = [
29             self._distance(self.v1, self.v2)^2,
30             self._distance(self.v2, self.v3)^2,
31             self._distance(self.v3, self.v1)^2
32         ]
33         # Sort the squared lengths
34         sides.sort()
35         # Check if the Pythagorean theorem holds
36         return sides[0] + sides[1] == sides[2]