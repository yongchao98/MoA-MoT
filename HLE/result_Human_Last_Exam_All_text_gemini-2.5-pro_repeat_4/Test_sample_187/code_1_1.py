13     def area(self):
14         # Heron's formula
15         a = self._distance(self.v1, self.v2)
16         b = self._distance(self.v2, self.v3)
17         c = self._distance(self.v3, self.v1)
18         s = (a + b + c) / 2
19         return sqrt(s * (s - a) * (s - b) * (s - c))