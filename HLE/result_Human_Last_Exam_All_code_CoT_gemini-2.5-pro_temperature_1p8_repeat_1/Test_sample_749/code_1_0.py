import math

def solve_branching_walk():
    """
    Calculates the limit of the probability that site 0 is visited
    by infinitely many particles in a branching random walk.
    """
    # Step 1: Define the jump probabilities based on the site color.
    p_right_if_red = 1/5
    p_left_if_red = 4/5
    p_right_if_blue = 4/5
    p_left_if_blue = 1/5

    # Step 2: Analyze the system in the limit h -> 0.
    # The probability is P(site is red) = h.
    # The branching probability is h.
    # We are interested in the limit h -> 0. Let's set h=0 to find the limit behavior.
    h = 0

    # Step 3: Calculate the average jump probabilities in the limit h -> 0.
    # p_avg is the average probability of a particle jumping right.
    # p_avg(h) = h * p_right_if_red + (1-h) * p_right_if_blue
    p_avg = h * p_right_if_red + (1 - h) * p_right_if_blue
    
    # q_avg is the average probability of a particle jumping left.
    # q_avg(h) = h * p_left_if_red + (1-h) * p_left_if_blue
    q_avg = h * p_left_if_red + (1 - h) * p_left_if_blue

    # Step 4: Define and calculate the key quantity K(h) in the limit h -> 0.
    # The criterion for the particle population at a fixed site to be self-sustaining
    # is K(h) >= 1, where K(h) is given by:
    # K(h) = (E[number of offspring])^2 * 4 * p_avg(h) * q_avg(h)
    # The number of "offspring" of a particle is either 1 (with prob 1-h) or 2 (with prob h).
    # The particle itself always survives, so the number of particles descending from one is either 1 or 2.
    # Let's consider the growth rate of lineages.
    # Another approach uses the condition that E[number of particles at origin] should not decay.
    # This leads to the condition (1+h)^2 * 2 * sqrt(p_avg*q_avg) >= 1
    # Let's calculate the value of this expression for h=0.
    
    K_0 = (1 + h)**2 * 2 * math.sqrt(p_avg * q_avg)
    
    # Step 5: Print the argument and the calculation.
    print("We analyze the behavior of the particle system as h -> 0.")
    print("The key question is whether the cloud of particles drifts away from site 0 or not.")
    print("This can be determined by a quantity K(h), which compares the particle creation rate to their movement.")
    print("If K(h) < 1, the particle cloud drifts to infinity, and any site is visited by only a finite number of particles.")
    print("This would mean the probability of infinitely many visits is 0.\n")

    print("Let's calculate K(h) = (1+h)^2 * 2 * sqrt(p_avg(h) * q_avg(h)) in the limit h -> 0.")
    print("First, we find the average jump probabilities p_avg and q_avg as h -> 0.")
    
    # Print equation for p_avg
    print(f"p_avg(h) = h * ({p_right_if_red}) + (1-h) * ({p_right_if_blue})")
    # Print equation for q_avg
    print(f"q_avg(h) = h * ({p_left_if_red}) + (1-h) * ({p_left_if_blue})\n")

    print(f"As h approaches 0:")
    print(f"p_avg(0) = 0 * {p_right_if_red} + (1-0) * {p_right_if_blue} = {p_avg}")
    print(f"q_avg(0) = 0 * {p_left_if_red} + (1-0) * {p_left_if_blue} = {q_avg}\n")
    
    print("Now we compute K(0):")
    p0_val = 4/5
    q0_val = 1/5
    term_sqrt_val = math.sqrt(p0_val * q0_val)
    print(f"K(0) = (1 + 0)^2 * 2 * sqrt(p_avg(0) * q_avg(0))")
    print(f"K(0) = (1)^2 * 2 * sqrt({p0_val} * {q0_val})")
    print(f"K(0) = 1 * 2 * sqrt({p0_val * q0_val})")
    print(f"K(0) = 2 * {term_sqrt_val}")
    print(f"K(0) = {K_0}\n")
    
    print(f"Since K(0) = {K_0:.2f} < 1, the condition for drift is met for small h.")
    print("This implies the entire particle cloud has a positive velocity to the right and moves towards +infinity.")
    print("Therefore, site 0 is visited by only a finite number of particles almost surely.")
    print("The limit of the probability is 0.")

solve_branching_walk()
<<<0>>>