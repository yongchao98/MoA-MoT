from math import sqrt, cos, sin, atan2 as atan

class Vector:
    def __init__(self, cx, cy):
        self.compx = cx
        self.compy = cy
        self.mag = sqrt(abs(self.compx**2 + self.compy**2))
        self.angle = atan(self.compy, self.compx)

    def calcDist(self, res, alt):
        at = sqrt(abs(alt/4.9))
        ot = sqrt(abs(float(res)**2 - self.mag**2))
        if at < ot:
            return (self.mag / 2) * at
        else:
            return (self.mag / 2) * ot

    def genDistCoor(self, res, alt):
        vd = self.calcDist(res, alt)
        return [vd*cos(self.angle), vd*sin(self.angle)]

def main_solution(cx, cy, res, alt):
    vector = Vector(cx, cy)
    result = vector.genDistCoor(res, alt)
    return result

# Given input
input_data = {'cx': -83.07813863062816, 'cy': 30.10843090078771, 'res': 91.90239164644267, 'alt': 136.16991674093114}
output = main_solution(input_data['cx'], input_data['cy'], input_data['res'], input_data['alt'])
print(output)