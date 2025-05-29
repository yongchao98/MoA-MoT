import numpy as np

def find_feasible_input():
    # Set a random seed for reproducibility
    np.random.seed(42)
    
    # Adjust constraints and global best values
    constraints = [[0, 10], [-5, 5], [0, 15], [0, 15], [0, 15]]
    global_best = [8, -3, 11, 12, 10]
    B = 0.5
    a = 0.1
    
    # Create a particle with the given constraints
    class Particle():
        def __init__(self, constraints):
            self.constraints = constraints
            self.pts  = np.zeros(len(constraints), dtype="float")
            self.randomize()

        def randomize(self):
            for i, (lowerbound, upperbound) in enumerate(self.constraints):
                self.pts[i]  = np.random.uniform(lowerbound, upperbound)

        def APSO(self, global_best, B, a):
            for i, pt in enumerate(self.pts):
                mu, sigma = 0, 1
                e = np.random.normal(mu, sigma)
                c = self.constraints[i]
                L = abs(c[1]-c[0])
                self.pts[i] = (1-B)*L*pt + B*L*global_best[i] + a*L*e

    particle = Particle(constraints)
    particle.APSO(global_best, B, a)
    
    return particle.pts.tolist()

print(find_feasible_input())