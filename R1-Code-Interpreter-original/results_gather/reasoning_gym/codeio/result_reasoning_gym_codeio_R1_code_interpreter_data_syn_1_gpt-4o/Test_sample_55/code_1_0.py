import math
import random

class Vector2:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def dot(self, other):
        return self.x * other.x + self.y * other.y

class Noise:
    def __init__(self):
        self.permutation = self.makePermutation()

    def shuffle(self, arrayToShuffle):
        for e in range(len(arrayToShuffle) - 1, 0, -1):
            index = math.floor(random.random() * (e + 1))
            temp = arrayToShuffle[e]
            arrayToShuffle[e] = arrayToShuffle[index]
            arrayToShuffle[index] = temp

    def makePermutation(self):
        permutation = []
        for i in range(256):
            permutation.append(i)
        self.shuffle(permutation)
        for i in range(256):
            permutation.append(permutation[i])
        return permutation

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

def find_coordinates(target_value, tolerance=0.01, search_range=20, step=0.1):
    noise_generator = Noise()
    x = 0
    while x < search_range:
        y = 0
        while y < search_range:
            noise_value = noise_generator.noise2D(x, y)
            if abs(noise_value - target_value) < tolerance:
                return x, y
            y += step
        x += step
    return None, None

x, y = find_coordinates(0.06191176711152624)
print(x, y)