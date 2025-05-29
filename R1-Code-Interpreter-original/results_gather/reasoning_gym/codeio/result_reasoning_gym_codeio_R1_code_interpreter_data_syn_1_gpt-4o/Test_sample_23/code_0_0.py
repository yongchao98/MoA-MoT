import math

class Setup(object):
    def __init__(self, echelle, articulation, sol_position):
        self.articulations = []
        self.membres = []
        self.echelle = echelle
        self.articulation = articulation
        self.sol_position = sol_position

    def scale(self, x: float, y: float) -> (float, float):
        xx = x * self.echelle
        yy = (1 - (y + self.sol_position)) * self.echelle
        return xx, yy

class Point(object):
    def __init__(self, setup: Setup, x: float, y: float, color="orange"):
        self.setup = setup
        self.x = x
        self.y = y
        self.color = color

    def rotate(self, ref, angle: float):
        x = self.x - ref.x
        y = self.y - ref.y
        r = math.sqrt(x * x + y * y)
        try:
            _a = math.acos(x / r)
            if y < 0:
                _a = 2*math.pi - _a
            _a += angle
            x = r * math.cos(_a) + ref.x
            y = r * math.sin(_a) + ref.y
            self.x = x
            self.y = y
        except:
            pass

class Articulation(Point):
    def __init__(self, setup: Setup, x0: float, y0: float, color="orange"):
        super().__init__(setup=setup, x=x0, y=y0, color=color)
        setup.articulations.append(self)

class Membre(object):
    def __init__(self, setup: Setup, longueur: float, art1: Articulation, art2: Articulation, masse: float = 0.0, color="red"):
        self.longueur = float(longueur)
        self.masse = masse
        self.art1 = art1
        self.art2 = art2
        self.color = color
        self.setup = setup
        self.setup.membres.append(self)

    def check_longueur(self) -> bool:
        longueur = math.sqrt((self.art1.x - self.art2.x)*(self.art1.x - self.art2.x) +
                             (self.art1.y - self.art2.y)*(self.art1.y - self.art2.y))
        return abs((longueur - self.longueur)/self.longueur) < 0.0001

class Body(object):
    def __init__(self, setup: Setup):
        self.setup = setup
        longueur_tibia = 0.25
        longueur_femur = 0.25
        longueur_tronc = 0.35

        self.tete = Articulation(setup=setup, x0=0.5, y0=longueur_tibia + longueur_femur + longueur_tronc)
        self.hanche = Articulation(setup=setup, x0=0.5, y0=longueur_tibia + longueur_femur)
        self.genou1 = Articulation(setup=setup, x0=0.5, y0=longueur_tibia)
        self.genou2 = Articulation(setup=setup, x0=0.5, y0=longueur_tibia, color="green")
        self.cheville1 = Articulation(setup=setup, x0=0.5, y0=0)
        self.cheville2 = Articulation(setup=setup, x0=0.5, y0=0)

        self.tronc = Membre(setup=setup, longueur=longueur_tronc, art1=self.tete, art2=self.hanche, masse=1)
        self.femur1 = Membre(setup=setup, longueur=longueur_femur, art1=self.hanche, art2=self.genou1, masse=1)
        self.tibia1 = Membre(setup=setup, longueur=longueur_tibia, art1=self.genou1, art2=self.cheville1, masse=1)
        self.femur2 = Membre(setup=setup, longueur=longueur_femur, art1=self.hanche, art2=self.genou2, masse=1, color="blue")
        self.tibia2 = Membre(setup=setup, longueur=longueur_tibia, art1=self.genou2, art2=self.cheville2, masse=1, color="blue")

        self.genou1.y = self.tibia1.longueur
        self.hanche.y = self.tibia1.longueur + self.femur1.longueur
        self.genou2.y = self.tibia2.longueur

def main_solution(echelle: float, articulation: float, sol_position: float, angle: float):
    setup = Setup(echelle, articulation, sol_position)
    body = Body(setup)

    # Rotate the body by the given angle
    for articulation in setup.articulations:
        articulation.rotate(body.hanche, angle)

    # Calculate the final positions of the articulations
    final_positions = [(art.x, art.y) for art in setup.articulations]

    return final_positions

# Given input
echelle = 1.9818248035457011
articulation = 0.17849871262899353
sol_position = 0.49888877376444907
angle = 1.9018452712119887

# Calculate the output
output = main_solution(echelle, articulation, sol_position, angle)
print(output)