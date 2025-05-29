import math

class Vector2:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def dot(self, other):
        return self.x * other.x + self.y * other.y

class Noise:
    def __init__(self):
        self.permutation = self.makePermutation()

    def makePermutation(self):
        # Use a fixed permutation table for determinism
        permutation = [151, 160, 137, 91, 90, 15,
                       131, 13, 201, 95, 96, 53, 194, 233, 7, 225,
                       140, 36, 103, 30, 69, 142, 8, 99, 37, 240,
                       21, 10, 23, 190, 6, 148, 247, 120, 234, 75,
                       0, 26, 197, 62, 94, 252, 219, 203, 117, 35,
                       11, 32, 57, 177, 33, 88, 237, 149, 56, 87,
                       174, 20, 125, 136, 171, 168, 68, 175, 74, 165,
                       71, 134, 139, 48, 27, 166, 77, 146, 158, 231,
                       83, 111, 229, 122, 60, 211, 133, 230, 220, 105,
                       92, 41, 55, 46, 245, 40, 244, 102, 143, 54,
                       65, 25, 63, 161, 1, 216, 80, 73, 209, 76,
                       132, 187, 208, 89, 18, 169, 200, 196, 135, 130,
                       116, 188, 159, 86, 164, 100, 109, 198, 173, 186,
                       3, 64, 52, 217, 226, 250, 124, 123, 5, 202,
                       38, 147, 118, 126, 255, 82, 85, 212, 207, 206,
                       59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
                       223, 183, 170, 213, 119, 248, 152, 2, 44, 154,
                       163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
                       129, 22, 39, 253, 19, 98, 108, 110, 79, 113,
                       224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
                       251, 34, 242, 193, 238, 210, 144, 12, 191, 179,
                       162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
                       49, 192, 214, 31, 181, 199, 106, 157, 184, 84,
                       204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
                       138, 236, 205, 93, 222, 114, 67, 29, 24, 72,
                       243, 141, 128, 195, 78, 66, 215, 61, 156, 180]
        return permutation + permutation

    def getConstantVector(self, v):
        h = v & 3
        if h == 0:
            return Vector2(1.0, 1.0)
        elif h == 1:
            return Vector2(-1.0, 1.0)
        elif h == 2:
            return Vector2(-1.0, -1.0)
        else:
            return Vector2(1.0, -1.0)

    def fade(self, t):
        return ((6 * t - 15) * t + 10) * t * t * t

    def lerp(self, t, a1, a2):
        return a1 + t * (a2 - a1)

    def noise2D(self, x, y):
        X = math.floor(x) & 255
        Y = math.floor(y) & 255
        xf = x - math.floor(x)
        yf = y - math.floor(y)
        topRight = Vector2(xf - 1.0, yf - 1.0)
        topLeft = Vector2(xf, yf - 1.0)
        bottomRight = Vector2(xf - 1.0, yf)
        bottomLeft = Vector2(xf, yf)
        valueTopRight = self.permutation[self.permutation[X + 1] + Y + 1]
        valueTopLeft = self.permutation[self.permutation[X] + Y + 1]
        valueBottomRight = self.permutation[self.permutation[X + 1] + Y]
        valueBottomLeft = self.permutation[self.permutation[X] + Y]
        dotTopRight = topRight.dot(self.getConstantVector(valueTopRight))
        dotTopLeft = topLeft.dot(self.getConstantVector(valueTopLeft))
        dotBottomRight = bottomRight.dot(self.getConstantVector(valueBottomRight))
        dotBottomLeft = bottomLeft.dot(self.getConstantVector(valueBottomLeft))
        u = self.fade(xf)
        v = self.fade(yf)
        return self.lerp(u, self.lerp(v, dotBottomLeft, dotTopLeft), self.lerp(v, dotBottomRight, dotTopRight))

def find_coordinates(target_value, step=0.001, limit=100):
    noise_generator = Noise()
    for x in range(-limit, limit):
        for y in range(-limit, limit):
            x_float = x * step
            y_float = y * step
            noise_value = noise_generator.noise2D(x_float, y_float)
            if math.isclose(noise_value, target_value, abs_tol=1e-6):
                return x_float, y_float
    return None, None

x, y = find_coordinates(0.06191176711152624)
print(x, y)