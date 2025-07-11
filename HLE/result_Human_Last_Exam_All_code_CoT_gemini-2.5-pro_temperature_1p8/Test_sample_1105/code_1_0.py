import math
import scipy.constants

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    given an initial velocity and a constant applied force.
    """
    # Physical constants and initial conditions
    # These are example values. You can change them to fit a specific scenario.
    G = scipy.constants.G  # Gravitational constant
    m = 1.0e12  # mass of asteroid in kg
    M = 5.0e4   # mass of spaceship in kg
    l0 = 10000.0 # initial distance in meters (10 km)
    v0 = 0.05    # initial velocity in m/s
    F = 0.01     # additional applied force in Newtons

    # Print initial parameters
    print("Given Parameters:")
    print(f"  Asteroid Mass (m): {m:.1e} kg")
    print(f"  Spaceship Mass (M): {M:.1e} kg")
    print(f"  Initial Distance (l0): {l0:.1f} m")
    print(f"  Initial Velocity (v0): {v0:.2f} m/s")
    print(f"  Applied Force (F): {F:.2f} N")
    print("-" * 30)

    # Derived term for convenience
    GMm = G * M * m

    # Check for escape conditions where l_max would be infinite
    # 1. If applied force F is greater than or equal to gravitational pull at l0
    force_gravity_at_l0 = GMm / (l0**2)
    if F >= force_gravity_at_l0:
        print("The applied force F is strong enough to overcome gravity from the start.")
        print("The spaceship will escape to infinity.")
        return

    # 2. If the total energy is enough to overcome the potential energy barrier
    # The peak of the potential energy barrier is at U_max = -2*sqrt(F*G*M*m)
    U_max = -2 * math.sqrt(F * GMm)
    # The initial energy is E_initial = 1/2*M*v0^2 - (G*M*m/l0 + F*l0)
    E_initial = 0.5 * M * v0**2 - (GMm / l0 + F * l0)

    if E_initial > U_max:
        print("The initial kinetic energy is too high.")
        print("The spaceship will overcome the potential barrier and escape to infinity.")
        return

    # If conditions for a finite distance are met, solve the quadratic equation.
    # F * l_max^2 - C * l_max + GMm = 0
    # where C = F*l0 + GMm/l0 - 0.5*M*v0^2

    C = F * l0 + (GMm / l0) - 0.5 * M * v0**2
    discriminant = C**2 - 4 * F * GMm

    if discriminant < 0:
         # This case is already covered by the energy check, but included for completeness.
        print("Mathematical solution is complex. This implies the spaceship escapes to infinity.")
        return

    # The physically correct solution is the smaller root
    l_max = (C - math.sqrt(discriminant)) / (2 * F)

    # Print the equation and the final answer
    print("The maximum distance l_max is the solution to the quadratic equation:")
    print(f"({F:.4f}) * l_max^2 - ({C:.2f}) * l_max + ({GMm:.2f}) = 0")
    print("-" * 30)
    print(f"The calculated maximum distance is: {l_max:.2f} meters.")


if __name__ == '__main__':
    calculate_max_distance()
    l_max = ( (0.01 * 10000.0 + (6.6743e-11 * 5.0e4 * 1.0e12 / 10000.0) - 0.5 * 5.0e4 * 0.05**2) - math.sqrt( ( (0.01 * 10000.0 + (6.6743e-11 * 5.0e4 * 1.0e12 / 10000.0) - 0.5 * 5.0e4 * 0.05**2) )**2 - 4 * 0.01 * (6.6743e-11 * 5.0e4 * 1.0e12)) ) / (2 * 0.01)
    final_answer = ( ( (0.01*10000.0) + (333715000000.0*50000.0/10000.0) - 0.5 * 50000.0 * 0.05**2 ) - math.sqrt( ( (0.01*10000.0) + (333715000000.0*50000.0/10000.0) - 0.5 * 50000.0 * 0.05**2 )**2 - 4*0.01* (333715000000.0*50000.0)  ) ) / (2 * 0.01)
    lmax_value = (371.21 - math.sqrt(371.21**2 - 4*0.01*3337150.0)) / (2*0.01)
    final_answer_value = 15277.53
    <<<l_max = ( (F*l0 + (G*M*m / l0) - 0.5*M*v0**2) - sqrt( (F*l0 + (G*M*m / l0) - 0.5*M*v0**2)**2 - 4*F*G*M*m) ) / (2*F)>>>